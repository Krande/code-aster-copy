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

from code_aster.Cata.Language.DataStructure import fonction_sdaster
from code_aster.Cata.Language.Syntax import _F, FACT, MACRO, SIMP
from code_aster.Commands import DEFI_CONSTANTE, DEFI_LIST_REEL
from code_aster.Supervis import UserMacro
from code_aster import CA

CA.init("--test", ERREUR=_F(ERREUR_F="EXCEPTION", ALARME="EXCEPTION"))
test = CA.TestCase()

# Check for consistency between the 'max' attribute and the returned value (#34412)


def check_ops(self, A1=None, A2=None, O1=None, O2=None, F1=None, F2=None, INFO=2):
    """Command mockup to check keywords conversion according to the 'max' attribute."""

    def verb(msg):
        if INFO > 1:
            print(msg)

    test.assertIsNone(self._result)
    if A1 is not None:
        verb(f"A1={A1}")
        test.assertIsInstance(A1, float, msg="A1 should be a single float")
        test.assertNotIn(type(A1), (list, tuple), msg="A1 should not be a list or tuple")
    if A2 is not None:
        verb(f"A2={A2}")
        test.assertNotIsInstance(A2, float, msg="A2 should not be a single float")
        test.assertIn(type(A2), (list, tuple), msg="A2 should be a list or tuple")
    if O1 is not None:
        verb(f"O1={O1}")
        test.assertIsInstance(O1, CA.Function, msg="O1 should be a single Function")
        test.assertNotIn(type(O1), (list, tuple), msg="O1 should not be a list or tuple")
    if O2 is not None:
        verb(f"O2={O2}")
        test.assertNotIsInstance(O2, CA.Function, msg="O2 should not be a single Function")
        test.assertIn(type(O2), (list, tuple), msg="O2 should be a list or tuple")
    if F1 is not None:
        verb(f"F1={F1}")
        test.assertIsInstance(F1, dict, msg="F1 should be a single dict")
        test.assertNotIn(type(F1), (list, tuple), msg="F1 should not be a list or tuple")
    if F2 is not None:
        verb(f"F2={F2}")
        test.assertNotIsInstance(F2, dict, msg="F2 should not be a single dict")
        test.assertIn(type(F2), (list, tuple), msg="F2 should be a list or tuple")

    return None


CHECK_cata = MACRO(
    nom="CHECK",
    sd_prod=None,
    A1=SIMP(statut="f", typ="R", max=1),
    A2=SIMP(statut="f", typ="R", max="**"),
    O1=SIMP(statut="f", typ=fonction_sdaster, max=1),
    O2=SIMP(statut="f", typ=fonction_sdaster, max=2),
    F1=FACT(statut="f", max=1, S1=SIMP(statut="o", typ="R", max=1)),
    F2=FACT(statut="f", max="**", S2=SIMP(statut="o", typ="R", max=1)),
)

CHECK = UserMacro("CHECK", CHECK_cata, check_ops)

cst = DEFI_CONSTANTE(VALE=1.0)

CHECK(A1=1.0)
CHECK(A1=(1.0,))
with test.assertRaisesRegex(CA.AsterError, "(At most 1 value|SUPERVIS_99)"):
    CHECK(A1=(1.0, 2.0))

CHECK(A2=1.0)
CHECK(A2=(1.0,))
CHECK(A2=(1.0, 2.0))

CHECK(O1=cst)
CHECK(O1=(cst,))
with test.assertRaisesRegex(CA.AsterError, "(At most 1 value|SUPERVIS_99)"):
    CHECK(O1=(cst, cst))

CHECK(O2=cst)
CHECK(O2=(cst,))
CHECK(O2=(cst, cst))

CHECK(F1=_F(S1=1.0))
CHECK(F1={"S1": 1.0})
CHECK(F1=(_F(S1=1.0),))
CHECK(F1=[{"S1": 1.0}])
with test.assertRaisesRegex(CA.AsterError, "(at most 1 occurrence|SUPERVIS_99)"):
    CHECK(F1=(_F(S1=1.0), _F(S1=2.0)))

CHECK(F2=_F(S2=1.0))
CHECK(F2={"S2": 1.0})
CHECK(F2=(_F(S2=1.0),))
CHECK(F2=[{"S2": 1.0}])
CHECK(F2=(_F(S2=1.0), _F(S2=2.0)))

# Check for float to int conversion (#34499)
nbpas = 4
nbarch = nbpas / 2

# accept float if int(value) == value
list0 = DEFI_LIST_REEL(DEBUT=0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=nbarch))
values = list0.getValues()

test.assertEqual(len(values), round(nbarch) + 1, msg="check length")

with test.assertRaisesRegex(CA.AsterError, "(NOMBRE.*Unexpected type|SUPERVIS_99)"):
    DEFI_LIST_REEL(DEBUT=0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=2.000001))

test.printSummary()

CA.close()
