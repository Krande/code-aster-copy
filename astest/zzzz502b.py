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

# Creation du maillage
monMaillage = CA.Mesh()

# Relecture du fichier MED
monMaillage.readMedFile("zzzz503a.mmed")

# Definition du modele Aster
monModel = CA.Model(monMaillage)
monModel.addModelingOnMesh(CA.Physics.Mechanics, CA.Modelings.Tridimensional)

monModel.build()

charCine = CA.MechanicalDirichletBC(monModel)
charCine.addBCOnCells(CA.PhysicalQuantityComponent.Dx, 0.0, "Bas")
test.assertEqual(charCine.getType(), "CHAR_CINE_MECA")

# Impossible d'affecter un blocage en temperature sur un DEPL
with test.assertRaises(RuntimeError):
    charCine.addBCOnCells(CA.PhysicalQuantityComponent.Temp, 0.0, "Haut")

charCine.build()

CA.close()
