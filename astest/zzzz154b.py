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


from code_aster.Cata.Syntax import AU_MOINS_UN, EXCLUS, FACT, MACRO, PRESENT_PRESENT, SIMP, UN_PARMI
from code_aster.Commands import *
from code_aster import CA
from code_aster.Supervis.ExecuteCommand import UserMacro

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))


class Spy:
    """Class to spy the keywords actually used."""

    keywords = {}

    @classmethod
    def sorted(cls):
        return sorted(cls.keywords)


def macro_test_ops(self, **kwargs):
    """Macro de test."""
    Spy.keywords = kwargs
    print("keywords are:", kwargs)


MACRO_TEST_CATA = MACRO(
    nom="MACRO_TEST",
    op=macro_test_ops,
    regles=(
        PRESENT_PRESENT("B", "A"),
        EXCLUS("A", "C"),
        AU_MOINS_UN("B", "C", "D", C=3.0, D=4.0),
        # EXCLUS("A", "C"),
        # rules order is important: MACRO_TEST(A=1) is not valid after setting the default
    ),
    A=SIMP(statut="f", typ="R"),
    B=SIMP(statut="f", typ="R"),
    C=SIMP(statut="f", typ="R"),
    D=SIMP(statut="f", typ="R"),
    AFFE=FACT(
        regles=(UN_PARMI("I", "J", I=0),),
        statut="f",
        max="**",
        NOM=SIMP(statut="o", typ="TXM"),
        I=SIMP(statut="f", typ="I"),
        J=SIMP(statut="f", typ="I"),
    ),
)

MACRO_TEST = UserMacro("MACRO_TEST", MACRO_TEST_CATA, macro_test_ops)


test = CA.TestCase()

MACRO_TEST(A=1.0)
test.assertSequenceEqual(["A", "C", "D"], Spy.sorted())

with test.assertRaises(CA.AsterError):
    MACRO_TEST(B=2.0)

MACRO_TEST(A=1.0, B=2.0)
test.assertSequenceEqual(["A", "B"], Spy.sorted())

with test.assertRaises(CA.AsterError):
    MACRO_TEST(A=1.0, C=3.0)

MACRO_TEST()
test.assertSequenceEqual(["C", "D"], Spy.sorted())

MACRO_TEST(AFFE=(_F(NOM="string")))
test.assertSequenceEqual(["AFFE", "C", "D"], Spy.sorted())
test.assertSequenceEqual(["I", "NOM"], sorted(Spy.keywords["AFFE"][0]))

MACRO_TEST(AFFE=(_F(NOM="string", J=10)))
test.assertSequenceEqual(["AFFE", "C", "D"], Spy.sorted())
test.assertSequenceEqual(["J", "NOM"], sorted(Spy.keywords["AFFE"][0]))

with test.assertRaises(CA.AsterError):
    MACRO_TEST(AFFE=(_F(NOM="string", I=0, J=10)))

test.printSummary()

CA.close()
