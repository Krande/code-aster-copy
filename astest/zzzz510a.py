# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

import numpy as N
import code_aster
from code_aster.Commands import *
from code_aster import MPI

code_aster.init("--test")

test = code_aster.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()

graph = code_aster.CommGraph()
balancer = code_aster.ObjectBalancer()
a = [i + rank * 10 for i in range(10)]
print(rank, a)
if rank == 0:
    graph.addCommunication(1)
    graph.addCommunication(3)
    balancer.addElementarySend(1, [0, 2])
    balancer.addElementarySend(3, [1, 3])
    balancer.setElementsToKeep([3])
elif rank == 1:
    graph.addCommunication(0)
    graph.addCommunication(2)
    balancer.addElementarySend(0, [5, 6])
    balancer.addElementarySend(2, [2, 4, 7])
elif rank == 2:
    graph.addCommunication(3)
    balancer.addElementarySend(3, [1, 7])
elif rank == 3:
    graph.addCommunication(1)
    balancer.addElementarySend(1, [8])

balancer.endElementarySendDefinition()
balancer.prepareCommunications()
result = balancer.balanceVectorOverProcesses(a)
print("Result ", result)

if rank == 0:
    test.assertEqual(result[0], 3.0)
    test.assertEqual(result[1], 4.0)
    test.assertEqual(result[2], 5.0)
    test.assertEqual(result[3], 6.0)
    test.assertEqual(result[4], 7.0)
    test.assertEqual(result[5], 8.0)
    test.assertEqual(result[6], 9.0)
    test.assertEqual(result[7], 15.0)
    test.assertEqual(result[8], 16.0)
elif rank == 1:
    test.assertEqual(result[0], 10.0)
    test.assertEqual(result[1], 11.0)
    test.assertEqual(result[2], 13.0)
    test.assertEqual(result[3], 18.0)
    test.assertEqual(result[4], 19.0)
    test.assertEqual(result[5], 0.0)
    test.assertEqual(result[6], 2.0)
    test.assertEqual(result[7], 38.0)
elif rank == 2:
    test.assertEqual(result[0], 20.0)
    test.assertEqual(result[1], 22.0)
    test.assertEqual(result[2], 23.0)
    test.assertEqual(result[3], 24.0)
    test.assertEqual(result[4], 25.0)
    test.assertEqual(result[5], 26.0)
    test.assertEqual(result[6], 28.0)
    test.assertEqual(result[7], 29.0)
    test.assertEqual(result[8], 12.0)
    test.assertEqual(result[9], 14.0)
    test.assertEqual(result[10], 17.0)
elif rank == 3:
    test.assertEqual(result[0], 30.0)
    test.assertEqual(result[1], 31.0)
    test.assertEqual(result[2], 32.0)
    test.assertEqual(result[3], 33.0)
    test.assertEqual(result[4], 34.0)
    test.assertEqual(result[5], 35.0)
    test.assertEqual(result[6], 36.0)
    test.assertEqual(result[7], 37.0)
    test.assertEqual(result[8], 39.0)
    test.assertEqual(result[9], 21.0)
    test.assertEqual(result[10], 27.0)
    test.assertEqual(result[11], 1.0)
    test.assertEqual(result[12], 3.0)

bMesh = code_aster.MeshBalancer()
if rank == 0:
    myMesh = code_aster.Mesh()
    myMesh.readMedFile("fort.20")
    bMesh.buildFromBaseMesh(myMesh)
    outMesh = bMesh.applyBalancingStrategy([1, 2, 9, 11, 17, 19, 25, 31])
elif rank == 1:
    outMesh = bMesh.applyBalancingStrategy([5, 6, 13, 15, 18, 20, 26, 32])
elif rank == 2:
    outMesh = bMesh.applyBalancingStrategy([7, 8, 14, 16, 22, 24, 28, 30])
elif rank == 3:
    outMesh = bMesh.applyBalancingStrategy([3, 4, 10, 12, 21, 23, 27, 29])

import medcoupling as medc, numpy as np

size = MPI.ASTER_COMM_WORLD.Get_size()
ttps = medc.GetUMeshGlobalInfo("fort.20", "Maillage_1")[0]
nodeNb = medc.GetUMeshGlobalInfo("fort.20", "Maillage_1")[3]
nbNodesProc = int(nodeNb / size)
startNode = nbNodesProc * rank
endNode = nbNodesProc * (rank + 1)
if rank == size - 1:
    endNode = nodeNb

coords = medc.MEDFileUMesh.LoadPartCoords(
    "fort.20", "Maillage_1", -1, -1, ["X", "Y", "Z"], startNode, endNode
)


nbNodes = endNode - startNode
nodes_mc = coords[0].toNumPyArray()
nodes_mc.shape = (nbNodes, 3)
coords = np.zeros(shape=(nbNodes, 3))
coords[:, :3] = nodes_mc[:, :3]

params = []
cts = []
for tps in ttps:
    for cellType, nbCellsType in tps:
        slc = medc.DataArray.GetSlice(slice(0, nbCellsType, 1), rank, size)
        params += [slc.start, slc.stop, slc.step]
        cts.append(cellType)
        pass
mrs = medc.MEDFileMeshReadSelector()
mrs.setNumberOfCoordsLoadSessions(10)
medFileUMesh = medc.MEDFileUMesh.LoadPartOf("fort.20", "Maillage_1", cts, params, -1, -1, mrs)

from code_aster.ObjectsExt.medctoaster import (
    ConvertToMEDFileGeoType,
    MEDCouplingMeshHelper,
    toAsterGeoType,
)

mesh = medFileUMesh
non_empty_levels = mesh.getNonEmptyLevels()
level_shift = 0  # Variable pour la creation d'une numérotation globale

# Sortie groupes
# Groupes de cells
groups_of_cells = {}
groups_of_nodes = {}

connectivity_aster = []
types = []

# Loop sur les niveaux
for level in sorted(non_empty_levels):
    mesh_level = mesh[level]
    # bis = mesh_level.getGroups()
    groups = mesh.getGroupsOnSpecifiedLev(level)
    for group in groups:
        if len(group) > 24:
            UTMESS("A", "MED_7")
            continue
        group_cells = mesh.getGroupArr(level, group)
        group_cells += 1
        group_cells += level_shift
        groups_of_cells.setdefault(group, [])
        groups_of_cells[group].extend(group_cells.toNumPyArray())

    # Loop sur toutes les types du niveau
    types_at_level = mesh_level.getAllGeoTypesSorted()
    for medcoupling_cell_type in types_at_level:

        # Cells du meme type
        cells_current_type = mesh_level.giveCellsWithType(medcoupling_cell_type)

        # Maillage 1 Single Geo Type (connectivité simple)
        mesh_current_type = medc.MEDCoupling1SGTUMesh(mesh_level[cells_current_type])

        # Type MED
        med_current_type = ConvertToMEDFileGeoType(medcoupling_cell_type)

        # Nombre de noeuds du type
        number_of_nodes_current_type = medc.MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(
            medcoupling_cell_type
        )

        # Remise en forme en tableau
        connectivity_current_type = (
            mesh_current_type.getNodalConnectivity()
            .toNumPyArray()
            .reshape(mesh_current_type.getNumberOfCells(), number_of_nodes_current_type)
        )

        if medcoupling_cell_type != medc.NORM_POINT1:
            connectivity_current_type = connectivity_current_type[
                :, MEDCouplingMeshHelper.getConnectivityMedToAster(medcoupling_cell_type)
            ]

        # Shift de 1
        connectivity_current_type += 1

        # Sauvegarde de la connectivité aster au km
        connectivity_aster.extend(connectivity_current_type)
        types.extend([med_current_type] * mesh_current_type.getNumberOfCells())
    level_shift += mesh_level.getNumberOfCells()

groups = mesh.getGroupsOnSpecifiedLev(1)
for group in groups:
    if len(group) > 24:
        UTMESS("A", "MED_7")
        continue
    group_nodes = mesh.getGroupArr(1, group).deepCopy()
    if not group_nodes:
        continue
    group_nodes += 1
    groups_of_nodes.setdefault(group, [])
    groups_of_nodes[group].extend(group_nodes.toNumPyArray())

# print("connectivity_aster ", connectivity_aster)

connect = [cell.tolist() for cell in connectivity_aster]
typ = [toAsterGeoType(i) for i in types]
coo = coords.flatten().tolist()

from code_aster import Mesh

meshOut = Mesh()
meshOut._initDefinition(3, coo, connect, typ, len(groups_of_cells), len(groups_of_nodes))

groups = list(zip(*groups_of_cells.items()))
if groups:
    meshOut._addGroupsOfCells(*groups)
groups = list(zip(*groups_of_nodes.items()))
if groups:
    meshOut._addGroupsOfNodes(*groups)

meshOut.debugPrint()

mesh2 = code_aster.IncompleteMesh()
mesh2.readMedFile("fort.20")
mesh2.debugPrint()

FIN()
