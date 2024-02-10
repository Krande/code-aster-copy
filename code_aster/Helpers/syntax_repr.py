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
from ..Supervis.visitors import GenericVisitor
from ..Utilities import force_list

REQ = chr(9670)
OPT = chr(9671)
DEF = chr(10192)
XOR = "/"
ALT = "|"
EMPTY = chr(8709)
INDENT = "    "


def _symbol(status):
    return {"o": REQ, "f": OPT, "d": DEF}.get(status)


def _typ2name(typ):
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


class BaseLine:
    """Base object"""

    def __init__(self, lvl, name) -> None:
        self.lvl = lvl
        self.name = name
        self._hide = None

    @property
    def hidden(self):
        return self._hide


class KwdLine(BaseLine):
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
        self.attrs = []  # no attrs means hidden
        self.into = []
        self.typ = None
        self.default = None
        self.cols = []  # maybe column for attrs, name, into/var

    @property
    def hidden(self):
        self._hide = not bool(self.attrs) or self.name == "reuse"
        return super().hidden

    def set(self, defs):
        symb = _symbol(defs.get("statut"))
        if not symb:
            return
        self.attrs.append(symb)
        self.typ = defs.get("typ")
        self.into = defs.get("into", [])
        self.default = defs.get("defaut")

    def _repr_value(self, value):
        if self.typ == "TXM":
            return f'"{value}"'
        return str(value)

    def repr(self):
        prefix = INDENT * self.lvl + " ".join(self.attrs) + " " + self.name + " = "
        if not self.into:
            try:
                value = _var(self.typ)
            except AttributeError:
                raise TypeError(self.name, self.typ)
            if self.default:
                value += f" (défaut: {self.default})"
            value = [self._repr_value(value)]
        else:
            values = list(self.into)
            # if not self.default and len(self.into) == 1:
            #     values.append(EMPTY)
            value = [f"{XOR} {self._repr_value(i)}" for i in values]
            value = []
            for i in self.into:
                value.append(f"{XOR} {self._repr_value(i)}")
                if self.default is not None and i == self.default:
                    value[-1] += f" (par défaut)"
            if not self.default and len(self.into) == 1:
                value[-1] += " (ou non renseigné)"
        lines = [prefix + value.pop(0) + ","]
        for remain in value:
            lines.append(" " * len(prefix) + remain + ",")
        return os.linesep.join(lines)


class FactLine(BaseLine):
    """For factor keyword:

    .. code-block::text

        AFFE = _F(
            ...
    """

    def repr(self):
        return INDENT * self.lvl + self.name + " = _F("


class CmdLine(BaseLine):
    """For command:

    .. code-block::text

        maillage = COMMAND(
            ...
    """

    def repr(self):
        return INDENT * self.lvl + self.name + "("


class CloseLine(BaseLine):
    """End of a block."""

    def __init__(self, lvl, end=",") -> None:
        super().__init__(lvl, None)
        self.end = end

    def repr(self):
        return INDENT * self.lvl + ")" + self.end


class DocSyntaxVisitor(GenericVisitor):
    """Visitor to the syntax of a Command"""

    def __init__(self, command):
        super().__init__()
        self.lines = []
        self.command = command
        self.mcsimp = None
        self.mcfact = None
        self.indent = 0

    def _visitComposite(self, step, userDict=None):
        """Visit a composite object (containing BLOC, FACT and SIMP objects)"""
        # loop first on keywords in rules
        # store "rules attrs"...?
        for name, entity in step.entities.items():
            if entity.getCataTypeId() == IDS.simp:
                self.mcsimp = name
            elif entity.getCataTypeId() == IDS.fact:
                self.mcfact = name
            elif entity.getCataTypeId() == IDS.bloc:
                self.bloc = name
            entity.accept(self)

    def visitCommand(self, step, userDict=None):
        """Visit a Command object"""
        # + reuse + sd_prod
        line = CmdLine(self.indent, self.command)
        self._append(line)
        self.indent += 1
        self._visitComposite(step, userDict)
        self.indent -= 1
        self._append(CloseLine(self.indent, end=""))

    def visitBloc(self, step, userDict=None):
        """Visit a Bloc object"""
        print(self.bloc)
        self._visitComposite(step, userDict)

    def visitFactorKeyword(self, step, userDict=None):
        """Visit a FactorKeyword object"""
        line = FactLine(self.indent, self.mcfact)
        self._append(line)
        self.indent += 1
        self._visitComposite(step, userDict)
        self.indent -= 1
        self._append(CloseLine(self.indent))

    def visitSimpleKeyword(self, step, skwValue):
        """Visit a SimpleKeyword object"""
        line = KwdLine(self.indent, self.mcsimp)
        line.set(step.definition)
        self._append(line)
        # print(line.repr())

    def _append(self, line):
        self.lines.append(line)

    def repr(self):
        """Text representation"""
        return os.linesep.join([i.repr() for i in self.lines if not i.hidden])


def testing(cmd):
    syntax = DocSyntaxVisitor(cmd.name)
    cmd.accept(syntax)
    print(syntax.repr())
