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

CA.init("--test")
test = CA.TestCase()

# Creation du maillage
mesh = CA.Mesh()
test.assertEqual(mesh.getType(), "MAILLAGE_SDASTER")

# Relecture du fichier MED
mesh.readMedFile("zzzz501a.mmed")

# help(mesh)

coord = mesh.getCoordinates()
test.assertEqual(coord.getType(), "CHAM_GEOM")
# help(coord)

coord1 = coord.copy()
test.assertEqual(coord.getValues(), coord1.getValues())
# check readonly access
print("coord[3] ", coord[3])
test.assertSequenceEqual(coord[3], [0.06666666666666667, 0.0, 0.0])
node = coord.getNode(3)
test.assertSequenceEqual(coord[3], [node.x(), node.y(), node.z()])

# Definition du modele Aster
model = CA.Model(mesh)
test.assertEqual(model.getType(), "MODELE_SDASTER")
model.addModelingOnMesh(CA.Physics.Mechanics, CA.Modelings.Tridimensional)

model.setSplittingMethod(CA.ModelSplitingMethod.GroupOfCells)
test.assertEqual(model.getSplittingMethod(), CA.ModelSplitingMethod.GroupOfCells)

model.setSplittingMethod(CA.ModelSplitingMethod.Centralized)
test.assertEqual(model.getSplittingMethod(), CA.ModelSplitingMethod.Centralized)

model.build()

# Definition du modele Aster
model2 = CA.Model(mesh)

with test.assertRaisesRegex(RuntimeError, "not allowed"):
    model2.addModelingOnMesh(CA.Physics.Thermal, CA.Modelings.DKT)

# Verification du comptage de référence sur le maillage
del mesh

with test.assertRaises(NameError):
    mesh

mesh2 = model.getMesh()
test.assertTrue("Tout" in mesh2.getGroupsOfCells())

# Vérification du debug
mesh2.debugPrint(66)

del coord

CA.saveObjects()

test.printSummary()

CA.close()
