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

from code_aster.Commands import *
from code_aster import CA

CA.init("--test", "--continue", ERREUR=_F(ALARME="EXCEPTION"))
test = CA.TestCase()

# test for usage of dict-like objects
test.assertEqual(len(ther_dict), 2, msg="check len()")
test.assertSequenceEqual(sorted(ther_dict.keys()), ["12", "13"], msg="check keys()")
ch1 = ther_dict["12"].getField("TEMP", 1)
# ch1.debugPrint()
ch3 = ther_dict["13"].getField("TEMP", 2)
# ch3.debugPrint()


NOM_CAS = ["12", "13"]
for cas in NOM_CAS:
    resu_cas = LIRE_RESU(
        TYPE_RESU="EVOL_THER",
        FORMAT="MED",
        MAILLAGE=mesh,
        TOUT_ORDRE="OUI",
        FORMAT_MED=_F(NOM_CHAM="TEMP", NOM_RESU=cas),
        UNITE=80,
    )
    test.assertEqual(resu_cas.getNumberOfIndexes(), 2)


test.assertAlmostEqual(max(ch1.getValues()), 1.0, msg="check ther1")
test.assertAlmostEqual(max(ch3.getValues()), 3.0, msg="check ther3")

result_ab = EXTR_CONCEPT(DICT=dict_test, NOM="ab")
result_ac = EXTR_CONCEPT(DICT=dict_test, NOM="ac")
cha = result_ab.getField("TEMP", 1)
chc = result_ac.getField("TEMP", 2)
test.assertAlmostEqual(max(cha.getValues()), 1.0, msg="check ther_ab")
test.assertAlmostEqual(max(chc.getValues()), 3.0, msg="check ther_ac")

test.printSummary()

CA.close()
