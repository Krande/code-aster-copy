# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
This module represents Commands syntax for the documentation.
"""

import os

from ..Cata.Language.SyntaxObjects import IDS
from ..Cata.Language import Rules
from ..Cata.Language import DataStructure as LDS
from ..Supervis.visitors import GenericVisitor
from ..Utilities import force_list

REQ = chr(9670)
OPT = chr(9671)
DEF = chr(10192)
XOR = "/"
ALT = "|"
# EMPTY = chr(8709)
INDENT = "    "


def _symbol(status):
    return {"o": REQ, "f": OPT, "d": DEF}.get(status)


def _typ2name(typ):
    if issubclass(typ, LDS.UnitBaseType):
        return "unit"
    name = typ.__name__.lower().replace("_sdaster", "")
    # exception for ParallelMesh
    if name.endswith("_p"):
        name = name.replace("_p", "")
    return name


def _var(typ):
    name = {"I": "int", "R": "float", "C": "complex", "TXM": "text"}.get(typ)
    if not name:
        names = set([_typ2name(i) for i in force_list(typ)])
        name = f" {XOR} ".join(names)
    return name


class Rule:
    def __init__(self, attrs, args) -> None:
        self._attrs = attrs
        self._args = args

    @classmethod
    def factory(cls, rule_object):
        if isinstance(rule_object, Rules.AtLeastOne):
            return Rule([REQ, ALT], rule_object.ruleArgs)
        if isinstance(rule_object, Rules.ExactlyOne):
            return Rule([REQ, XOR], rule_object.ruleArgs)
        if isinstance(rule_object, Rules.AtMostOne):
            return Rule([OPT, XOR], rule_object.ruleArgs)
        if isinstance(rule_object, Rules.OnlyFirstPresent):
            return Rule([XOR, XOR], rule_object.ruleArgs)
        raise NotImplementedError(rule_object)

    @property
    def args(self):
        return self._args

    def consumed(self):
        """Mark the rule as consumed/started"""
        assert self._attrs
        self._attrs[0] = " "

    def involved(self, keyword):
        """Tell if a keyword is involved in the rule"""
        return keyword in self._args


class BaseLine:
    """Base object"""

    def __init__(self, lvl, name) -> None:
        self._lvl = lvl
        self._name = name
        self._hide = None

    @property
    def hidden(self):
        return self._hide

    @property
    def offset(self):
        return INDENT * self._lvl

    def set(self, defs, rules):
        pass


class KwdLine(BaseLine):
    """Common for simple and factor keyworeds."""

    def __init__(self, lvl, name) -> None:
        super().__init__(lvl, name)
        self._attrs = []
        self._rules = []

    def set(self, defs, rules):
        symb = _symbol(defs.get("statut"))
        if not symb:
            return
        self._attrs.append(symb)
        self._rules = [rule for rule in rules if rule.involved(self._name)]

    @property
    def attrs(self):
        if self._rules:
            assert len(self._rules) == 1, f"{self._name}: {self._rules}"
            attrs = self._rules[0]._attrs[:]
            self._rules[0].consumed()
            return attrs
        return self._attrs

    @property
    def offset(self):
        return INDENT * self._lvl + " ".join(self.attrs) + " "


class SimpKwdLine(KwdLine):
    """Content for :

    .. code-block:: text

        ♦     FICHIER = / 'GLOBALE',
        attrs  name   =   into         type
        ♦     FICHIER = text,

    Replace "[I]", "[TXM]" by "int", "text" as "value". No type when into exists.
    + min/max, "l_int" ?
    """

    def __init__(self, lvl, name) -> None:
        super().__init__(lvl, name)
        self._into = []
        self._typ = None
        self._default = None
        self._cols = []  # maybe column for attrs, name, into/var

    @property
    def hidden(self):
        self._hide = not bool(self._attrs) or self._name == "reuse"
        return super().hidden

    def set(self, defs, rules):
        super().set(defs, rules)
        self._typ = defs.get("typ")
        self._into = defs.get("into", [])
        self._default = defs.get("defaut")

    def _repr_value(self, value):
        if self._typ == "TXM":
            return f'"{value}"'
        return str(value)

    def repr(self):
        prefix = self.offset + self._name + " = "
        if not self._into:
            try:
                value = _var(self._typ)
            except AttributeError:
                raise TypeError(self._name, self._typ)
            if self._default:
                value += f" (défaut: {self._default})"
            value = [self._repr_value(value)]
        else:
            values = list(self._into)
            value = [f"{XOR} {self._repr_value(i)}" for i in values]
            value = []
            for i in self._into:
                value.append(f"{XOR} {self._repr_value(i)}")
                if self._default is not None and i == self._default:
                    value[-1] += f" (par défaut)"
            if not self._default and len(self._into) == 1:
                value[-1] += " (ou non renseigné)"
        lines = [prefix + value.pop(0) + ","]
        for remain in value:
            lines.append(" " * len(prefix) + remain + ",")
        return os.linesep.join(lines)


class FactLine(KwdLine):
    """For factor keyword:

    .. code-block::text

        AFFE = _F(
            ...
    """

    def repr(self):
        return self.offset + self._name + " = _F("


class CmdLine(BaseLine):
    """For command:

    .. code-block::text

        maillage = COMMAND(
            ...
    """

    def repr(self):
        return self.offset + self._name + "("


class CondLine(BaseLine):
    """For conditional blocks"""

    def repr(self):
        return self.offset + "# " + self._name


class CloseLine(BaseLine):
    """End of a block."""

    def __init__(self, parent, end=",") -> None:
        super().__init__(0, "")
        self._parent = parent
        if isinstance(parent, CondLine):
            self._mark = ""
        else:
            self._mark = ")" + end

    @property
    def offset(self):
        return " " * len(self._parent.offset)

    def repr(self):
        return self.offset + self._mark


class DocSyntaxVisitor(GenericVisitor):
    """Visitor to the syntax of a Command"""

    def __init__(self, command):
        super().__init__()
        self._lines = []
        self._command = command
        self._mcsimp = None
        self._mcfact = None
        self._indent = 0
        self._rstack = []

    def _visitComposite(self, step, userDict=None):
        """Visit a composite object (containing BLOC, FACT and SIMP objects)"""
        extracted = step.entities
        entities = []
        rules = self._rstack[-1][:]
        while extracted:
            name, entity = extracted.popitem(last=False)
            entities.append((name, entity))
            for rule in rules:
                if rule.involved(name):
                    for kwd in rule.args:
                        ent = extracted.pop(kwd, None)
                        if ent:
                            entities.append((kwd, ent))

        for name, entity in entities:
            if entity.getCataTypeId() == IDS.simp:
                self._mcsimp = name
            elif entity.getCataTypeId() == IDS.fact:
                self._mcfact = name
            elif entity.getCataTypeId() == IDS.bloc:
                self._bloc = entity.getCondition()
            entity.accept(self)

    def _visitCompositeWrap(self, step, name, line_class, end, userDict=None):
        """Visit a Command object"""
        self._rstack.append([Rule.factory(rule) for rule in step.rules])
        line = line_class(self._indent, name)
        line.set(step.definition, self._rstack[-1])
        self._append(line)
        self._indent += 1
        self._visitComposite(step, userDict)
        self._indent -= 1
        self._append(CloseLine(line, end=end))
        self._rstack.pop()

    def visitCommand(self, step, userDict=None):
        """Visit a Command object"""
        # + reuse + sd_prod
        self._visitCompositeWrap(step, self._command, CmdLine, "", userDict)

    def visitBloc(self, step, userDict=None):
        """Visit a Bloc object"""
        self._visitCompositeWrap(step, self._bloc, CondLine, None, userDict)

    def visitFactorKeyword(self, step, userDict=None):
        """Visit a FactorKeyword object"""
        self._visitCompositeWrap(step, self._mcfact, FactLine, ",", userDict)

    def visitSimpleKeyword(self, step, skwValue):
        """Visit a SimpleKeyword object"""
        line = SimpKwdLine(self._indent, self._mcsimp)
        line.set(step.definition, self._rstack[-1])
        self._append(line)
        # print(line.repr())

    def _append(self, line):
        self._lines.append(line)

    def repr(self, legend=False):
        """Text representation"""
        text = [i.repr() for i in self._lines if not i.hidden]
        if legend:
            text.extend(
                [
                    "",
                    f"{REQ} : obligatoire",
                    f"{OPT} : optionnel",
                    f"{DEF} : présent par défaut",
                    f"{XOR} : un parmi",
                    f"{ALT} : plusieurs choix possibles",
                ]
            )
        text.append("")
        return os.linesep.join(text)


def testing(cmd):
    syntax = DocSyntaxVisitor(cmd.name)
    cmd.accept(syntax)
    print(syntax.repr(legend=True))
