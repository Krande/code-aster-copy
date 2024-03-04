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

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

TCPU = INFO_EXEC_ASTER(LISTE_INFO="TEMPS_RESTANT")
remaining = TCPU["TEMPS_RESTANT", 1]
test.assertGreater(remaining, 30.0)
test.assertLess(remaining, 9999.0)

tab = INFO_EXEC_ASTER(LISTE_INFO="ETAT_UNITE", FICHIER="file1")
state = tab["ETAT_UNITE", 1].strip()
test.assertEqual(state, "FERME")

unit = DEFI_FICHIER(ACTION="ASSOCIER", FICHIER="file1", ACCES="NEW")
test.assertEqual(unit, 99)

tab = INFO_EXEC_ASTER(LISTE_INFO="ETAT_UNITE", FICHIER="file1")
state = tab["ETAT_UNITE", 1].strip()
test.assertEqual(state, "OUVERT")

DEFI_FICHIER(ACTION="LIBERER", UNITE=unit)

tab = INFO_EXEC_ASTER(LISTE_INFO="ETAT_UNITE", FICHIER="file1")
state = tab["ETAT_UNITE", 1].strip()
test.assertEqual(state, "FERME")

test.printSummary()

FIN()
