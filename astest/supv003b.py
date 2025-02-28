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
from code_aster.Cata.Language.Syntax import FACT, MACRO, SIMP
from code_aster.Commands import DEFI_CONSTANTE
from code_aster.Supervis import UserMacro
from code_aster import CA

CA.init("--test", ERREUR=_F(ERREUR_F="EXCEPTION", ALARME="EXCEPTION"))
test = CA.TestCase()


def check_ops(self, A1=None, A2=None, O1=None, O2=None, INFO=2):
    """Command mockup to check keywords conversion according to the 'max' attribute."""

    def verb(msg):
        if INFO > 1:
            print(msg)

    test.assertIsNone(self._result)
    if A1 is not None:
        verb(f"{A1=}")
        test.assertIs(type(A1), float, msg="A1 should be a single float")
        test.assertNotIn(type(A1), (list, tuple), msg="A1 should not be a list or tuple")
    if A2 is not None:
        verb(f"{A2=}")
        test.assertIsNot(type(A2), float, msg="A2 should not be a single float")
        test.assertIn(type(A2), (list, tuple), msg="A2 should be a list or tuple")
    if O1 is not None:
        verb(f"{O1=}")
        test.assertIs(type(O1), CA.Function, msg="O1 should be a Function")
        test.assertNotIn(type(O1), (list, tuple), msg="O1 should not be a list or tuple")
    if O2 is not None:
        verb(f"{O2=}")
        test.assertIsNot(type(O2), CA.Function, msg="O2 should not be a single Function")
        test.assertIn(type(O2), (list, tuple), msg="O2 should be a list or tuple")

    return None


CHECK_cata = MACRO(
    nom="CHECK",
    sd_prod=None,
    A1=SIMP(statut="f", typ="R", max=1),
    A2=SIMP(statut="f", typ="R", max="**"),
    O1=SIMP(statut="f", typ=fonction_sdaster, max=1),
    O2=SIMP(statut="f", typ=fonction_sdaster, max=2),
)

CHECK = UserMacro("CHECK", CHECK_cata, check_ops)

cst = DEFI_CONSTANTE(VALE=1.0)

CHECK(A1=1.0)
CHECK(A1=(1.0,))
with test.assertRaisesRegex(CA.AsterError, "At most 1 value"):
    CHECK(A1=(1.0, 2.0))

CHECK(A2=1.0)
CHECK(A2=(1.0,))
CHECK(A2=(1.0, 2.0))

CHECK(O1=cst)
CHECK(O1=(cst,))
with test.assertRaisesRegex(CA.AsterError, "At most 1 value"):
    CHECK(O1=(cst, cst))

CHECK(O2=cst)
CHECK(O2=(cst,))
CHECK(O2=(cst, cst))

test.printSummary()

CA.close()
