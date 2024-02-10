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
import unittest
from pathlib import Path

from ..Cata import Commands as CMD
from ..Cata.Language import DataStructure as LDS
from ..Cata.Language import Rules
from ..Cata.Language.Syntax import FACT, PRESENT_ABSENT, SIMP
from ..Cata.Language.SyntaxObjects import IDS
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
    """Store Rules parameters."""

    def __init__(self, attrs, args, next=[]) -> None:
        self._attrs = attrs
        self._args = args
        self._first = True
        self._next = next

    @classmethod
    def factory(cls, rule_object):
        """Create an instance from a Cata Rule object."""
        if isinstance(rule_object, Rules.AtLeastOne):
            return Rule([REQ, ALT], rule_object.ruleArgs)
        if isinstance(rule_object, Rules.ExactlyOne):
            return Rule([REQ, XOR], rule_object.ruleArgs)
        if isinstance(rule_object, Rules.AtMostOne):
            return Rule([OPT, XOR], rule_object.ruleArgs)
        if isinstance(rule_object, Rules.OnlyFirstPresent):
            next = [ALT] if len(rule_object.ruleArgs) > 2 else []
            return Rule([OPT, XOR], rule_object.ruleArgs, next=next)
        raise NotImplementedError(rule_object)

    @property
    def args(self):
        """Return the involved keywords."""
        return self._args

    def consumed(self):
        """Mark the rule as consumed/started"""
        assert self._attrs
        for i, symb in enumerate(self._attrs[:-1]):
            if symb.strip():
                self._attrs[i] = " "
                break
        if self._first:
            self._attrs.extend(self._next)
        self._first = False

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
        """Tell if the block should be hidden"""
        return self._hide

    @property
    def offset(self):
        """Current offset/indentation"""
        return INDENT * self._lvl

    def set(self, defs, rules):
        """Set block settings from the keyword definition"""
        pass


class KwdLine(BaseLine):
    """Common for simple and factor keyworeds."""

    def __init__(self, lvl, name) -> None:
        super().__init__(lvl, name)
        self._attrs = []
        self._rules = []

    def set(self, defs, rules):
        """Set block settings from the keyword definition"""
        symb = _symbol(defs.get("statut"))
        if not symb:
            return
        self._attrs.append(symb)
        self._rules = [rule for rule in rules if rule.involved(self._name)]

    @property
    def attrs(self):
        """Return the current symbols."""
        if self._rules:
            if len(self._rules) > 1:
                print(f"WARNING: several rules for {self._name}")
            attrs = self._rules[0]._attrs[:]
            self._rules[0].consumed()
            return attrs
        return self._attrs

    @property
    def offset(self):
        """Current offset/indentation"""
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
        """Tell if the block should be hidden"""
        self._hide = not bool(self._attrs) or self._name == "reuse"
        return super().hidden

    def set(self, defs, rules):
        """Set block settings from the keyword definition"""
        super().set(defs, rules)
        self._typ = defs.get("typ")
        self._into = defs.get("into", [])
        self._default = defs.get("defaut")

    def _repr_value(self, value):
        if self._typ == "TXM":
            return f'"{value}"'
        return str(value)

    def repr(self):
        """Representation of the block."""
        prefix = self.offset + self._name + " = "
        if not self._into:
            try:
                value = _var(self._typ)
            except AttributeError:
                raise TypeError(self._name, self._typ)
            if self._default:
                value += f" (défaut: {self._repr_value(self._default)})"
            value = [value]
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
        """Representation of the block."""
        return self.offset + self._name + " = _F("


class CmdLine(BaseLine):
    """For command:

    .. code-block::text

        maillage = COMMAND(
            ...
    """

    def repr(self):
        """Representation of the block."""
        return self.offset + self._name + "("


class CondLine(BaseLine):
    """For conditional blocks"""

    def repr(self):
        """Representation of the block."""
        return self.offset + "# " + self._name


class CloseLine(BaseLine):
    """End of a block."""

    def __init__(self, parent, end=",") -> None:
        super().__init__(0, "")
        self._parent = parent
        self._mark = ")" + end

    @property
    def hidden(self):
        """Tell if the block should be hidden"""
        return isinstance(self._parent, CondLine)

    @property
    def offset(self):
        """Current offset/indentation"""
        return " " * len(self._parent.offset)

    def repr(self):
        """Representation of the block."""
        return self.offset + self._mark


class DocSyntaxVisitor:
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

    def visitMacro(self, step, userDict=None):
        """Visit a MacroCommand object"""
        self.visitCommand(step, userDict)

    def visitBloc(self, step, userDict=None):
        """Visit a Bloc object"""
        self._visitCompositeWrap(step, self._bloc, CondLine, "", userDict)

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

    def repr(self, legend=False, codeblock=False):
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
        text = os.linesep.join(text)
        if codeblock:
            text = [INDENT + line for line in text.splitlines()]
            text.insert(0, "")
            text.insert(0, ".. code-block:: text")
            text = os.linesep.join(text)
        return text


class TestDoc(unittest.TestCase):
    """Test for DocSyntaxVisitor"""

    output = None
    # uncomment the following line to write outputs of DocSyntax
    output = Path("/tmp/syntax")

    @classmethod
    def _testcmd(cls, cmd):
        syntax = DocSyntaxVisitor(cmd.name)
        cmd.accept(syntax)
        text = syntax.repr(legend=True, codeblock=True)
        if cls.output:
            cls.output.mkdir(parents=True, exist_ok=True)
            with open(cls.output / (cmd.name.lower() + ".rst"), "w") as frst:
                frst.write(text)
        return text

    def test00_lire_maillage(self):
        text = self._testcmd(CMD.LIRE_MAILLAGE)

    def test00_affe_materiau(self):
        self._testcmd(CMD.AFFE_MATERIAU)

    def test00_calc_champ(self):
        self._testcmd(CMD.CALC_CHAMP)

    def test01_onlyfirst(self):
        fact = FACT(
            statut="d",
            regles=(
                PRESENT_ABSENT(
                    "RESI_REFE_RELA", "RESI_GLOB_MAXI", "RESI_GLOB_RELA", "RESI_COMP_RELA"
                ),
            ),
            RESI_REFE_RELA=SIMP(statut="f", typ="R"),
            RESI_GLOB_MAXI=SIMP(statut="f", typ="R"),
            RESI_GLOB_RELA=SIMP(statut="f", typ="R"),
            RESI_COMP_RELA=SIMP(statut="f", typ="R"),
            ITER_GLOB_MAXI=SIMP(statut="f", typ="I", defaut=10),
            ITER_GLOB_ELAS=SIMP(statut="f", typ="I", defaut=25),
        )
        syntax = DocSyntaxVisitor("TEST")
        syntax._mcfact = "CONVERGENCE"
        fact.accept(syntax)
        text = syntax.repr(legend=False)
        # print(text)
        self.assertTrue(f"{DEF} CONVERGENCE = _F(" in text, msg="CONVERGENCE")
        self.assertTrue(f"{OPT} / RESI_REFE_RELA" in text, msg="RESI_REFE_RELA nook")
        self.assertTrue("/ | RESI_GLOB_MAXI" in text, msg="RESI_GLOB_MAXI nook")
        self.assertTrue("  | RESI_GLOB_RELA" in text, msg="RESI_GLOB_RELA nook")
        self.assertTrue("  | RESI_COMP_RELA" in text, msg="RESI_COMP_RELA nook")
