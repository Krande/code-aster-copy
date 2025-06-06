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

test = CA.TestCase()

try:
    DEBUT(CODE="OUI", ERREUR=_F(ALARME="EXCEPTION", ERREUR_F="EXCEPTION"))

except CA.AsterError as exc:
    assert False, "no exception should be thrown."

print("Checking translation of AsterErrorCpp as AsterError(Py) and catched as AsterError...")
try:
    CA.raiseAsterError("SUPERVIS_2")
    assert False, "This line should not be reached as an exception should have been thrown by now."
except CA.AsterError as exc:
    test.assertEqual(exc.id_message, "SUPERVIS_2")

with test.assertRaisesRegex(CA.AsterError, "commandes DEBUT et POURSUITE"):
    CA.raiseAsterError("SUPERVIS_2")


print("Checking translation of AsterErrorCpp as AsterError(Py) and catched as Exception...")
try:
    CA.raiseAsterError("SUPERVIS_2")
    assert False, "This line should not be reached as an exception should have been thrown by now."
except Exception as exc:
    test.assertEqual(exc.id_message, "SUPERVIS_2")

with test.assertRaisesRegex(Exception, "commandes DEBUT et POURSUITE"):
    CA.raiseAsterError("SUPERVIS_2")


print("Checking AsterErrorCpp thrown from C++ methods and catched as AsterError...")
mesh = CA.Mesh()
try:
    mesh.readAsterFile("fort.18")
except CA.AsterError as exc:
    test.assertEqual(exc.id_message, "MODELISA_1")

print("Checking AsterErrorCpp thrown from Fortran operator and catched as AsterError...")
try:
    LIRE_MAILLAGE(UNITE=55, FORMAT="ASTER")
except CA.AsterError as exc:
    test.assertEqual(exc.id_message, "MODELISA4_93")


test.printSummary()

FIN()

FIN()
