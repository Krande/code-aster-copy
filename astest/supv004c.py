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


from code_aster.Commands import *
from code_aster import CA
from code_aster.Utilities import ExecutionParameter, Options

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))


class MyCalc:

    """Example where commands calls are encapsulated in a Python objet."""

    _list = _formula = None

    def __init__(self, nb_step, val_max, formula):
        """Store command result as attribute."""
        self._list = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=val_max, NOMBRE=nb_step))
        self._formula = formula

    def get_formula(self):
        """Direct as returned value."""
        return FORMULE(VALE=self._formula, NOM_PARA="X")

    def create_table(self):
        """In method."""
        tabx = CREA_TABLE(LISTE=_F(LISTE_R=self._list.getValuesAsArray(), PARA="X"))

        result = CALC_TABLE(
            TABLE=tabx, ACTION=_F(OPERATION="OPER", FORMULE=self.get_formula(), NOM_PARA="X2")
        )
        return result


test = CA.TestCase()

mycalc = MyCalc(10, 100.0, "X * 2")

# example to hide commands syntax
ExecutionParameter().disable(Options.ShowSyntax)

tab = mycalc.create_table()

# check for dependencies
deps = tab.getDependencies()
test.assertEqual(len(deps), 0, msg="tab dependencies")

# restore printing of commands syntax
ExecutionParameter().enable(Options.ShowSyntax)

test.assertEqual(tab.userName, "result", msg="check name")

FIN()
