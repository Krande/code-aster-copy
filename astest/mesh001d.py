# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

import code_aster

code_aster.init("--test")

test = code_aster.TestCase()
mesh = code_aster.Mesh()
mesh.readMedFile("zzzz503a.mmed")

mcmesh = mesh.createMedCouplingMesh()

test.assertEqual(mcmesh.getName(),
                 mesh.getName(), mesh.getName())


test.assertEqual(mcmesh.getNumberOfNodes(),
                 mesh.getNumberOfNodes(), mesh.getNumberOfNodes())

test.assertEqual(sum(mcmesh.getNumberOfCellsAtLevel(i) for i in mcmesh.getNonEmptyLevels()),
                     mesh.getNumberOfCells(), mesh.getNumberOfCells())

test.assertEqual(mcmesh.getCoords().getValues(),
                 mesh.getCoordinates().getValues())


test.assertEqual(set(mcmesh.getGroupsOnSpecifiedLev(1)),
                 set(mesh.getGroupsOfNodes()))

test.assertEqual(set(mcmesh.getGroupsOnSpecifiedLev(-1)),
                 set(mesh.getGroupsOfCells()))

aster_connect = mesh.getMedConnectivity()
connect_seg2_aster =  [tuple(i) for i in aster_connect if len(i)==2]
connect_seg2_mc = mcmesh[-2].convertNodalConnectivityToStaticGeoTypeMesh()
connect_seg2_mc.rearrange(2)
connect_seg2_mc+=1

test.assertEqual(connect_seg2_mc.getValuesAsTuple(),
                 connect_seg2_aster)

connect_quad4_aster =  [tuple(i) for i in aster_connect if len(i)==4]
connect_quad4_mc = mcmesh[-1].convertNodalConnectivityToStaticGeoTypeMesh()
connect_quad4_mc.rearrange(4)
connect_quad4_mc+=1
test.assertEqual(connect_quad4_mc.getValuesAsTuple(),
                 connect_quad4_aster)

connect_hexa8_aster =  [tuple(i) for i in aster_connect if len(i)==8]
connect_hexa8_mc = mcmesh[0].convertNodalConnectivityToStaticGeoTypeMesh()
connect_hexa8_mc.rearrange(8)
connect_hexa8_mc+=1
test.assertEqual(connect_hexa8_mc.getValuesAsTuple(),
                 connect_hexa8_aster)

test.printSummary()

code_aster.close()
