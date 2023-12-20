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
from code_aster.Cata.DataStructure import evol_ther_dict, evol_ther
from code_aster.Cata.Syntax import FACT, MACRO, SIMP
from code_aster.Supervis.ExecuteCommand import UserMacro


def macro_test_ops(self, AFFE, **kwargs):
    """Macro de test."""
    result = CA.ThermalResultDict()
    for fact in AFFE:
        result[fact["NOM_CAS"]] = fact["RESULTAT"]
    return result


MACRO_TEST_CATA = MACRO(
    nom="MACRO_TEST",
    op=macro_test_ops,
    sd_prod=evol_ther_dict,
    AFFE=FACT(
        statut="o",
        max="**",
        NOM_CAS=SIMP(statut="o", typ="TXM"),
        RESULTAT=SIMP(statut="o", typ=evol_ther),
    ),
)

MACRO_TEST = UserMacro("MACRO_TEST", MACRO_TEST_CATA, macro_test_ops)
