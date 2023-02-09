# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

import code_aster
from code_aster.Commands import *

code_aster.init("--test")

from zzzz154a_cmd import MACRO_TEST

test = code_aster.TestCase()

mesh = code_aster.Mesh.buildSquare()

ther1 = CREA_CHAMP(
    MAILLAGE=mesh,
    TYPE_CHAM="NOEU_TEMP_R",
    OPERATION="AFFE",
    AFFE=_F(NOM_CMP="TEMP", VALE=1.0, TOUT="OUI"),
)
ther2 = CREA_CHAMP(
    MAILLAGE=mesh,
    TYPE_CHAM="NOEU_TEMP_R",
    OPERATION="AFFE",
    AFFE=_F(NOM_CMP="TEMP", VALE=2.0, TOUT="OUI"),
)
ther3 = CREA_CHAMP(
    MAILLAGE=mesh,
    TYPE_CHAM="NOEU_TEMP_R",
    OPERATION="AFFE",
    AFFE=_F(NOM_CMP="TEMP", VALE=3.0, TOUT="OUI"),
)

result12 = CREA_RESU(
    OPERATION="AFFE",
    TYPE_RESU="EVOL_THER",
    NOM_CHAM="TEMP",
    AFFE=(_F(CHAM_GD=ther1, INST=1.0), _F(CHAM_GD=ther2, INST=2.0)),
)

result13 = CREA_RESU(
    OPERATION="AFFE",
    TYPE_RESU="EVOL_THER",
    NOM_CHAM="TEMP",
    AFFE=(_F(CHAM_GD=ther1, INST=1.0), _F(CHAM_GD=ther3, INST=3.0)),
)

ther_dict = code_aster.ThermalResultDict("ther_dict")
test.assertEqual(len(ther_dict), 0, msg="check len()")
test.assertEqual(ther_dict.getType(), "EVOL_THER_DICT", msg="check type")

ther_dict["12"] = result12
ther_dict["13"] = result13

test.assertEqual(len(ther_dict), 2, msg="check len()")
test.assertIsNone(ther_dict.get("unknown"), msg="check invalid access")
with test.assertRaises(KeyError):
    _ = ther_dict["unknown"]
with test.assertRaises(TypeError):
    ther_dict["invalid"] = ther1

# ch1 = ther_dict["12"].getField("TEMP", 1)
# ch3 = ther_dict["13"].getField("TEMP", 2)
# test.assertAlmostEqual(max(ch1.getValues()), 1.0, msg="check ther1")
# test.assertAlmostEqual(max(ch3.getValues()), 3.0, msg="check ther3")

# test for usage of dict-like objects in commands
dict_test = MACRO_TEST(
    AFFE=(_F(NOM_CAS="ab", RESULTAT=result12), _F(NOM_CAS="ac", RESULTAT=result13))
)

test.printSummary()

code_aster.close()
