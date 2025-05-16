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
Unittest for CommandSyntax.
"""


from code_aster.Commands import *
from code_aster import CA
from code_aster.Cata.Commands import DEBUT
from code_aster.Supervis import CommandSyntax

CA.init("--test", "--abort", "--debug", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

# just for testing import of command description from path
syntax = CommandSyntax("code_aster.Cata.Commands.DEBUT")
test.assertEqual(syntax.getName(), "DEBUT")

catadef = DEBUT.definition
test.assertIn("INFO", catadef.simple_keywords)
test.assertNotIn("INFO", catadef.factor_keywords)
test.assertNotIn("CODE", catadef.factor_keywords)
test.assertIn("CODE", catadef.simple_keywords)

userkw = _F(
    IMPR_MACRO="OUI",
    ERREUR=_F(ERREUR_F="EXCEPTION"),
    RESERVE_CPU=_F(POURCENTAGE=0.20),
    IGNORE_ALARM=("alarm1", "alarm22", "alarm333", "alarm4444"),
    INFO=2,
)

syntax.setResult("result", "type")
test.assertEqual(syntax.getResultName(), "result")
test.assertEqual(syntax.getResultType(), "type")
test.assertSequenceEqual(syntax.getres(), ("result", "type", "DEBUT"))


syntax.define(userkw)
test.assertEqual(syntax.getFactorKeywordNbOcc("ERREUR"), 1)

test.assertTrue(syntax.getexm("", "IMPR_MACRO"))
test.assertTrue(syntax.getexm(" ", "CODE"))
test.assertFalse(syntax.getexm("ERREUR", "ALARM"))
test.assertTrue(syntax.getexm("ERREUR", ""))
test.assertTrue(syntax.getexm("ERREUR", " "))

size, length = syntax.getltx("", "IGNORE_ALARM", 0, 3, 7)
test.assertEqual(size, -4)
test.assertSequenceEqual(length, [6, 7, 7])

size, length = syntax.getltx("", "LANG", 0, 3, 7)
test.assertEqual(size, 0)
test.assertSequenceEqual(length, [])

# fake test of getvid with strings
size, value = syntax.getvid("", "IGNORE_ALARM", 0, 3)
test.assertEqual(size, -4)
test.assertSequenceEqual(value, ["alarm1", "alarm22", "alarm333"])

size, value = syntax.getvtx("", "IGNORE_ALARM", 999, 6)
test.assertEqual(size, 4)
test.assertSequenceEqual(value, ["alarm1", "alarm22", "alarm333", "alarm4444"])

size, value = syntax.getvis("", "INFO", 999, 1)
test.assertEqual(size, 1)
test.assertSequenceEqual(value, [2])

size, value = syntax.getvr8("RESERVE_CPU", "POURCENTAGE", 0, 1)
test.assertEqual(size, 1)
test.assertSequenceEqual(value, [0.2])

kws, types = syntax.getmjm("RESERVE_CPU", 0, 99)
test.assertEqual(len(kws), len(types))
test.assertEqual(len(kws), 1)
test.assertSequenceEqual(kws, ["POURCENTAGE"])
test.assertSequenceEqual(types, ["R8"])

kws, types = syntax.getmjm(" ", 0, 99)
test.assertEqual(len(kws), len(types))
test.assertEqual(len(kws), 4)
test.assertSequenceEqual(kws, ["IGNORE_ALARM", "IMPR_MACRO", "INFO", "MODE"])
test.assertSequenceEqual(types, ["TX", "TX", "IS", "TX"])

syntax.free()

# some more complex cases (factor keywords under conditional blocks)
syntax = CommandSyntax("DYNA_VIBRA")
test.assertEqual(syntax.getName(), "DYNA_VIBRA")

test.assertTrue(syntax.getexm("", "TYPE_CALCUL"))
test.assertTrue(syntax.getexm("", "INFO"))

test.assertTrue(syntax.getexm("IMPRESSION", ""))
test.assertTrue(syntax.getexm("IMPRESSION", "NIVEAU"))

test.assertTrue(syntax.getexm("EXCIT", "CHARGE"))
test.assertTrue(syntax.getexm("EXCIT", "D_FONC_DT"))

syntax.free()

# case using a part of syntax
syntax = CommandSyntax("code_aster.Cata.Commons.c_comportement.C_COMPORTEMENT_MNL")
test.assertEqual(syntax.getName(), "C_COMPORTEMENT_MNL")

test.assertTrue(syntax.getexm("COMPORTEMENT", "TOUT"))
test.assertTrue(syntax.getexm("COMPORTEMENT", "POST_INCR"))

syntax.free()

test.printSummary()

CA.close(exit=True)

raise RuntimeError("no should pass here!")
