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
monMaillage = CA.Mesh()

# test de relecture d'un fichier Gmsh
monMaillage.readGmshFile("ssnv187a.msh")

# test du format Gibi
mtest = CA.Mesh()
mtest.readGibiFile("zzzz364a.mgib")

coord = monMaillage.getCoordinates()

# check readonly access
test.assertSequenceEqual(coord[3], [0.0, 1.0, 0.0])

# Definition du modele Aster
monModel = CA.Model(monMaillage)
monModel.addModelingOnMesh(CA.Physics.Mechanics, CA.Modelings.Tridimensional)

# delete monMaillage and check that the C++ object still exists because
# referenced by the model object
del monMaillage
with test.assertRaises(NameError):
    monMaillage.getCoordinates()

# delete/overwrite monModel, coord object still exists
monModel = 1
test.assertEqual(coord[3], [0.0, 1.0, 0.0])

test.printSummary()

CA.close()
