# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------

"""
This module defines the base objects that allow to use the legacy syntax
of code_aster commands.

Note:

    - If a keyword (simple or factor) is not provided by the user, its value
      is None.

    - If a factor keyword is present by default (``statut='d'``) but not filled by
      the user, its value is ``{}`` if ``max=1`` or ``[]`` if ``max > 1``.

    - Only one level of factor keywords is supported.

The keyword ``reuse`` is (will be) deprecated. So it must be present in catalogs
to import command files that contain it. But it is always optional (see
``checkMandatory``) because it is removed during import.
"""

import types
from collections import OrderedDict

from . import DataStructure as DS
from .DataStructure import DataStructure, UnitBaseType
from .SyntaxChecker import CheckerError, checkCommandSyntax
from .SyntaxUtils import (
    _F,
    add_none_sdprod,
    block_utils,
    debug_message2,
    disable_0key,
    enable_0key,
    force_list,
    is_undefined,
    sorted_dict,
    value_is_sequence,
)


class SyntaxId:
    """Container of the id of syntax objects.

    This list of type identifiers can be extended but never change between
    two releases of code_aster.
    """

    simp, fact, bloc, command = list(range(4))


IDS = SyntaxId()

UNDEF = object()


# Must stay identical to `AsterStudy.datamodel.general.ConversionLevel`
class ConversionLevel:
    """
    Enumerator for the level of conversion requirements.

    Attributes:
        NoFail: Do not fail, not *strict*.
        Naming: Requires that all command results are explicitly named.
        Type: Requires a valid type definition.
        Keyword: Requires that all keywords are valid.
        Syntaxic: Requires a valid syntax of all the commands.
        Restore: Requires a conversion without error during restore.
        Any: All conversion must pass.
        Partial: Allows to make a partial conversion (to be used with
            another level).
        NoGraphical: Force to load all stages in text mode.
    """

    NoFail = 0x00
    Naming = 0x01
    Type = 0x02
    Keyword = 0x04
    Syntaxic = Naming | Type | Keyword
    Restore = 0x08
    Any = Syntaxic | Restore
    Partial = 0x10
    NoGraphical = 0x20


class CataDefinition(OrderedDict):
    """Dictionary to store the definition of syntax objects.
    Iteration over the elements is ordered by type: SimpleKeyword, FactorKeyword
    and Bloc.
    """

    @property
    def entities(self):
        """Return the all entities.

        Returns:
            dict: dict of all entities (keywords and conditional blocks) of the
            object.
        """
        kws = self._filter_entities((SimpleKeyword, Bloc, FactorKeyword), with_block=False)
        return sorted_dict(kws)

    @property
    def keywords(self):
        """Return the simple and factor keywords contained in the object.

        The keywords are sorted in the order of definition in the catalog.
        This is a workaround waiting for Python 3.6 and integration
        of `PEP-0468`_.

        Returns:
            dict: dict of all simple and factor keywords of the object.

        .. _PEP-0468: https://www.python.org/dev/peps/pep-0468/
        """
        kws = self._filter_entities((SimpleKeyword, FactorKeyword))
        return sorted_dict(kws)

    @property
    def factor_keywords(self):
        """Return the factor keywords contained in the object.

        Returns:
            dict: dict of all factor keywords of the object.
        """
        kws = self._filter_entities((FactorKeyword,))
        return sorted_dict(kws)

    @property
    def simple_keywords(self):
        """Return the simple keywords contained in the object.

        Returns:
            dict: dict of all simple keywords of the object.
        """
        kws = self._filter_entities((SimpleKeyword,))
        return sorted_dict(kws)

    def _filter_entities(self, typeslist, with_block=True):
        """Filter entities by type recursively.

        Returns:
            dict: dict of entities of the requested types.
        """
        entities = CataDefinition()
        for key, value in list(self.items()):
            if isinstance(value, typeslist):
                entities[key] = value
            elif with_block and isinstance(value, Bloc):
                bloc_kwds = value.definition._filter_entities(typeslist)
                # update without overwriting factorkeywords but extending them
                for kwd, obj in bloc_kwds.items():
                    if type(obj) is FactorKeyword and kwd in entities:
                        # do not change existing FactorKeyword!
                        all_kwds = entities[kwd].definition.copy()
                        all_kwds.update(obj.definition)
                        entities[kwd] = FactorKeyword(all_kwds)
                    else:
                        entities[kwd] = obj
        return entities

    def iterItemsByType(self):
        """Iterator over dictionary's pairs with respecting order:
        SimpleKeyword, FactorKeyword and Bloc objects"""
        keysR = [k for k, v in list(self.items()) if type(v) is SimpleKeyword]
        keysI = [k for k, v in list(self.items()) if type(v) is FactorKeyword]
        keysS = [k for k, v in list(self.items()) if type(v) is Bloc]
        for key in keysR + keysI + keysS:
            yield (key, self[key])


class UIDMixing:
    """Sub class for UID based classes.

    Arguments:
        uid (int): Object's id.
    """

    _new_id = -1
    _id = None

    @classmethod
    def new_id(cls):
        UIDMixing._new_id += 1
        return UIDMixing._new_id

    def __init__(self):
        self._id = self.new_id()

    @property
    def uid(self):
        """Attribute that holds unique *id*"""
        return self._id

    def __lt__(self, other):
        if other is None or not hasattr(other, "uid"):
            return True
        return self._id < other.uid

    def __eq__(self, other):
        if other is None or not hasattr(other, "uid"):
            return False
        return self._id == other.uid


class PartOfSyntax(UIDMixing):
    """
    Generic object that describe a piece of syntax.
    """

    def __init__(self, curDict):
        """Initialization"""
        super().__init__()
        self._definition = CataDefinition(curDict)
        regles = curDict.get("regles")
        if regles and type(regles) not in (list, tuple):
            regles = (regles,)
        self._rules = regles or []

    def getCataTypeId(self):
        """Get the Cata type of object.
        Should be sub-classed.

        Returns:
            int: type id of Cata object: -1 if not defined.
        """
        return -1

    @property
    def name(self):
        """str: Name of the object."""
        return self._definition.get("nom", "")

    @property
    def udocstring(self):
        """unicode: Documentation of the object."""
        doc = self.docstring
        if type(doc) is not str:
            doc = str(doc, "utf-8", "replace")
        return doc

    @property
    def docstring(self):
        """str: Documentation of the object."""
        return self._definition.get("fr", "")

    @property
    def definition(self):
        """dict: Attribute containing the syntax definition"""
        return self._definition

    @property
    def rules(self):
        """dict: Attribute containing the list of rules"""
        return self._rules

    # : for backward compatibility
    regles = rules

    def __repr__(self):
        """Simple representation"""
        return "%s( %r )" % (self.__class__, self._definition)

    @property
    def entities(self):
        """Return the all entities contained in the object.

        Returns:
            dict: dict of all entities (keywords and conditional blocks) of the
            object.
        """
        return self._definition.entities

    # : for backward compatibility (and avoid changing `pre_seisme_nonl`)
    entites = entities

    @property
    def keywords(self):
        """Return the simple and factor keywords contained in the object.

        Returns:
            dict: dict of all simple and factor keywords of the object.
        """
        return self._definition.keywords

    def accept(self, visitor, syntax=None):
        """Called by a Visitor"""
        raise NotImplementedError("must be defined in a subclass")

    def _def_status(self):
        """Wrapper that returns the value of 'statut' after precondition."""
        definition = self.definition
        value = definition.get("statut", "f")
        if value == "c":
            return value
        # In AsterStudy, a simple keyword with a default value...
        if self.getCataTypeId() == IDS.simp and self.hasDefaultValue():
            try:
                # ... with derivated type of UnitBaseType...
                if issubclass(definition.get("typ"), UnitBaseType):
                    # ... is mandatory
                    return "o"
            except TypeError:
                pass
        return value

    def isMandatory(self):
        """Tell if this keyword is mandatory"""
        return self._def_status() == "o"

    def isOptional(self):
        """Tell if this keyword is optional"""
        return self._def_status() == "f"

    def isHidden(self):
        """Tell if this keyword should be hidden"""
        return self._def_status() == "c"

    def hasDefaultValue(self):
        """Tell if the keyword has a default value"""
        return False

    def addDefaultKeywords(self, userSyntax, _parent_ctxt=None):
        """Add default keywords into the user dict of keywords.

        Optional keywords that are not defined are set to None/undefined value.

        The values given in argument (userSyntax) preempt on the definition
        ones.

        Arguments:
            userSyntax (dict): dict of the keywords as filled by the user,
                **changed** in place.
            _parent_ctxt (dict): contains the keywords as known in the parent.
                This context is used to evaluate block conditions.
        """
        if userSyntax is None:
            return
        ctxt = _parent_ctxt.copy() if _parent_ctxt else {}
        for key, kwd in self.definition.iterItemsByType():
            if isinstance(kwd, SimpleKeyword):
                # use the default
                kwd.addDefaultKeywords(key, userSyntax, ctxt)
            elif isinstance(kwd, FactorKeyword):
                # if not given by the user, default value is None if optional
                # or {} if present by default
                userFact = userSyntax.get(key, kwd.defaultValue())
                if userFact is None:
                    pass
                elif type(userFact) in (list, tuple):
                    for userOcc in userFact:
                        try:
                            kwd.addDefaultKeywords(userOcc, ctxt)
                        except AttributeError:
                            raise CheckerError(
                                TypeError, f"{key}: Unexpected type: {type(userOcc)}", []
                            )
                else:
                    try:
                        kwd.addDefaultKeywords(userFact, ctxt)
                    except AttributeError:
                        raise CheckerError(
                            TypeError, f"{key}: Unexpected type: {type(userFact)}", []
                        )
                    if kwd.is_list():
                        userFact = [userFact]
                userSyntax[key] = userFact
                ctxt[key] = userSyntax[key]
            elif isinstance(kwd, Bloc):
                if kwd.isEnabled(ctxt):
                    kwd.addDefaultKeywords(userSyntax, ctxt)

    def checkMandatory(self, userSyntax, stack, _parent_ctxt=None):
        """Check that the mandatory keywords are provided by the user.

        Warning: Default keywords must be added before visiting the objects.
        Warning: this does not check recursively, only the current level.

        Arguments:
            userSyntax (dict): dict of the user and default keywords.
            stack (list): used to give contextual informations in error
                messages.
            _parent_ctxt (dict): contains the keywords as known in the parent.
                This context is used to evaluate block conditions.
        """
        ctxt = _parent_ctxt.copy() if _parent_ctxt else {}
        # update parent context with local keywords
        ctxt.update(userSyntax)

        for key, kwd in self.definition.iterItemsByType():
            if isinstance(kwd, (SimpleKeyword, FactorKeyword)):
                if key == "reuse":  # reuse is deprecated
                    continue
                    # raise KeyError("reuse has been removed")
                # pragma pylint: disable=no-member
                if kwd.isMandatory() and kwd.undefined(userSyntax.get(key)):
                    debug_message2("mandatory keyword =", key, ":", kwd)
                    debug_message2("given syntax =", userSyntax)
                    stack.append(key)
                    raise KeyError("Keyword {0} is mandatory".format(key))
            elif isinstance(kwd, Bloc):
                if kwd.isEnabled(ctxt):
                    kwd.checkMandatory(userSyntax, stack, ctxt)
            # else: sdprod, fr...

    def getKeyword(self, userKeyword, userSyntax, _parent_ctxt=None):
        """Return the keyword in the current composite object.

        Warning: Default keywords must be added before visiting the objects.

        Arguments:
            userKeyword (str): name of the searched keyword.
            userSyntax (dict): dict of the user and default keywords.
            _parent_ctxt (dict): contains the keywords as known in the parent.
                This context is used to evaluate block conditions.
        """
        # search in keywords list
        found = self.definition.get(userKeyword)
        if found:
            return found

        ctxt = _parent_ctxt.copy() if _parent_ctxt else {}
        # update parent context with local keywords
        ctxt.update(userSyntax)

        # search in BLOC objects
        for _, kwd in self.definition.iterItemsByType():
            if not isinstance(kwd, Bloc):
                continue
            # debug_message2("block", key, repr(kwd.getCondition()), "with", ctxt,
            #                "enabled ?", kwd.isEnabled(ctxt))
            if not kwd.isEnabled(ctxt):
                continue
            found = kwd.getKeyword(userKeyword, userSyntax, ctxt)
            if found:
                break
        return found

    def getRules(self, userSyntax, _parent_ctxt=None):
        """Return the rules to be applied to the given keywords.

        Warning: Default keywords must be added before visiting the objects.

        Arguments:
            userSyntax (dict): dict of the user and default keywords.
            _parent_ctxt (dict): contains the keywords as known in the parent.
                This context is used to evaluate block conditions.
        """
        ctxt = _parent_ctxt.copy() if _parent_ctxt else {}
        # update parent context with local keywords
        ctxt.update(userSyntax)

        rules = list(self.rules)
        # search in BLOC objects
        for _, kwd in self.definition.iterItemsByType():
            if not isinstance(kwd, Bloc):
                continue
            if not kwd.isEnabled(ctxt):
                continue
            rules.extend(kwd.getRules(userSyntax))
        return rules

    @classmethod
    def undefined(cls, value):
        """Return *True* if the value is a null value (undefined keyword),
        *False* otherwise."""
        return is_undefined(value)

    def is_list(self):
        """Tell if the value should be stored as list."""
        return self.definition.get("max", 1) != 1


class SimpleKeyword(PartOfSyntax):
    """
    Objet mot-clé simple équivalent à SIMP dans les capy
    """

    def __init__(self, curDict):
        """Initialization"""
        super(SimpleKeyword, self).__init__(curDict)
        if "val_min" in self._definition or "val_max" in self._definition:
            typ = self._definition["typ"]
            assert typ in ("I", "R", "C") or issubclass(
                typ, UnitBaseType
            ), "'val_min/val_max' not allowed for type" " '{0}'".format(typ)

    def getCataTypeId(self):
        """Get the type id of SimpleKeyword.

        Returns:
            int: type id of SimpleKeyword.
        """
        return IDS.simp

    def accept(self, visitor, syntax=None):
        """Called by a Visitor"""
        visitor.visitSimpleKeyword(self, syntax)

    def _context(self, value):
        """Print contextual informations"""
        debug_message2("CONTEXT: value={!r}, type={}".format(value, type(value)))
        debug_message2("CONTEXT: definition: {}".format(self))

    def hasDefaultValue(self):
        """Tell if the keyword has a default value"""
        return self.defaultValue(UNDEF) is not UNDEF

    def defaultValue(self, default=None):
        """Return the default value or *default*."""
        return self.definition.get("defaut", default)

    def addDefaultKeywords(self, key, userSyntax, _parent_ctxt):
        """Add the default value if not provided by the user dict.

        Arguments:
            userSyntax (dict): dict of the keywords as filled by the user,
                **changed** in place.
            _parent_ctxt (dict): contains the keywords as known in the parent.
        """
        value = userSyntax.get(key, self.defaultValue())
        # do not insert None value
        if value is None:
            return
        if self.is_list():
            if value is not None and not value_is_sequence(value):
                value = [value]
        userSyntax[key] = value
        _parent_ctxt[key] = userSyntax[key]


class FactorKeyword(PartOfSyntax):
    """
    Objet mot-clé facteur equivalent de FACT dans les capy
    """

    def getCataTypeId(self):
        """Get the type id of FactorKeyword.

        Returns:
            int: type id of FactorKeyword.
        """
        return IDS.fact

    def accept(self, visitor, syntax=None):
        """Called by a Visitor"""
        visitor.visitFactorKeyword(self, syntax)

    def defaultValue(self):
        """Return the *default* for this factor keyword.
        If the keyword is present by default, its default value is {} if max=1
        or [] if max > 1.
        If the keyword does not exist by default, its default value is None."""
        if self._def_status() not in ("c", "d"):
            return None
        if self.is_list():
            return [_F()]
        return _F()


class Bloc(PartOfSyntax):
    """
    Objet Bloc équivalent à BLOC dans les capy
    """

    def getCataTypeId(self):
        """Get the type id of Bloc.

        Returns:
            int: type id of Bloc.
        """
        return IDS.bloc

    def accept(self, visitor, syntax=None):
        """Called by a Visitor"""
        visitor.visitBloc(self, syntax)

    def getCondition(self):
        """Return the BLOC condition"""
        cond = self.definition.get("condition")
        assert cond is not None, "A bloc must have a condition!"
        return cond

    def isEnabled(self, context):
        """Tell if the block is enabled by the given context"""
        eval_context = {}
        eval_context.update(DS.__dict__)
        eval_context.update(block_utils(eval_context))
        # evaluate Python variables if present
        for key, value in list(context.items()):
            eval_context[key] = getattr(value, "evaluation", value)
        try:
            enabled = eval(self.getCondition(), {}, eval_context)
        except AssertionError:
            raise
        except Exception:
            # TODO: re-enable CataError, it seems me a catalog error!
            # raise CataError("Error evaluating {0!r}: {1}".format(
            #                 self.getCondition(), str(exc)))
            enabled = False
        return enabled


class Command(PartOfSyntax):
    """
    Object Command qui représente toute la syntaxe d'une commande
    """

    _call_callback = None

    def get_compat_syntax(self):
        """Return the `compat_syntax` keywords converter."""
        return self.definition.get("compat_syntax", lambda _: None)

    @classmethod
    def register_call_callback(cls, callback):
        """Register *callback* to be called in place of the default method.

        Register *None* to revert to default.

        Arguments:
            callback (callable): Function to call with signature:
                (*Command* instance, ``**kwargs``).
        """
        cls._call_callback = callback

    def getCataTypeId(self):
        """Get the type id of Command.

        Returns:
            int: type id of Command.
        """
        return IDS.command

    def can_reuse(self):
        """Tell if the result can be a reused one."""
        reentr = self.definition.get("reentrant", "").split(":")
        return reentr and reentr[0] in ("o", "f")

    def accept(self, visitor, syntax=None):
        """Called by a Visitor"""
        visitor.visitCommand(self, syntax)

    def call_default(self, **args):
        """Execute the command, only based on the command description.

        Keyword arguments:
            __strict__ (bool): If True, the syntax of the command must be valid.
                Otherwise the generic DataStructure type is returned.
        """
        strict = args.pop("__strict__", ConversionLevel.Syntaxic)
        if strict & ConversionLevel.Syntaxic:
            checkCommandSyntax(self, args)
            resultType = self.get_type_sd_prod(**args)
        else:
            try:
                resultType = self.get_type_sd_prod(**args)
            except:
                return DataStructure()

        if resultType is None:
            return None

        return resultType()

    def __call__(self, **args):
        """Simulate the command execution."""
        if Command._call_callback is None:
            return self.call_default(**args)
        else:
            return Command._call_callback(self, **args)

    def get_type_sd_prod(self, **ctxt):
        """Return the type of the command result."""
        resultType = self.definition.get("sd_prod")
        if resultType is not None:
            if type(resultType) is types.FunctionType:
                self.addDefaultKeywords(ctxt)
                # if os.environ.get('DEBUG'):
                # print "COMMAND:", self.name
                # print "CTX1:", ctxt.keys()
                # print "CTX1:", ctxt
                add_none_sdprod(resultType, ctxt)
                resultType = self.build_sd_prod(resultType, ctxt)
            return resultType

    def build_sd_prod(self, sdprodFunc, ctxt):
        """Call the `sd_prod` function"""
        enable_0key(ctxt)
        resultType = sdprodFunc(**ctxt)
        disable_0key(ctxt)
        return resultType

    def get_all_types(self):
        """Return the list of all possible types that the command can return.

        Returns:
            list[.DataStructure.DataStructure] or list[list[.DataStructure.DataStructure]]:
            List of possible
            types or a list of list if there are additional results for
            macro-commands.
        """
        sd_prod = self.definition.get("sd_prod")
        if type(sd_prod) is types.FunctionType:
            args = {}
            add_none_sdprod(sd_prod, args)
            args["__all__"] = True
            return force_list(sd_prod(**args))
        else:
            return (sd_prod,)


class Operator(Command):
    pass


class Macro(Command):
    """Specialization of the Command object for MACRO"""

    def build_sd_prod(self, sdprodFunc, ctxt):
        """Call the `sd_prod` function"""
        enable_0key(ctxt)
        if "self" in ctxt:
            del ctxt["self"]
        resultType = sdprodFunc(self, **ctxt)
        disable_0key(ctxt)
        return resultType

    def type_sdprod(self, result, astype):
        """Define the type of the result."""
        if result is None:
            return
        result.settype(astype)

    def accept(self, visitor, syntax=None):
        """Called by a Visitor"""
        visitor.visitMacro(self, syntax)


class Procedure(Command):
    pass


class Formule(Macro):
    pass


class CataError(Exception):
    """Exception raised in case of error in the catalog."""

    pass
