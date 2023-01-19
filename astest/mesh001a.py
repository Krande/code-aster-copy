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

import code_aster
from code_aster.Commands import RECU_TABLE
import numpy as np

code_aster.init("--test")

# check Mesh object API
test = code_aster.TestCase()

# from MED format
mesh = code_aster.Mesh()
mesh.readMedFile("zzzz503a.mmed", verbose=2)
with test.assertRaisesRegex(AssertionError, "not empty"):
    mesh.readMedFile("zzzz503a.mmed")

test.assertFalse(mesh.isParallel())
test.assertFalse(mesh.isQuadratic())
test.assertEqual(mesh.getDimension(), 3)
test.assertEqual(mesh.getNumberOfNodes(), 27)
test.assertEqual(mesh.getNumberOfCells(), 56)

# test dimension
TABG = RECU_TABLE(CO=mesh, NOM_TABLE="CARA_GEOM")
test.assertAlmostEqual(TABG["X_MIN", 1], 0.0)
test.assertAlmostEqual(TABG["X_MAX", 1], 1.0)
test.assertAlmostEqual(TABG["Y_MIN", 1], 0.0)
test.assertAlmostEqual(TABG["Y_MAX", 1], 1.0)
test.assertAlmostEqual(TABG["Z_MIN", 1], 0.0)
test.assertAlmostEqual(TABG["Z_MAX", 1], 1.0)
test.assertAlmostEqual(TABG["AR_MIN", 1], 0.5)
test.assertAlmostEqual(TABG["AR_MAX", 1], 0.5)

# test groups
# do the same thing (compatibily with ParallelMesh)
test.assertSequenceEqual(
    sorted(mesh.getGroupsOfNodes()), ["A", "B", "Bas", "C", "D", "E", "F", "G", "H", "Haut"]
)
test.assertSequenceEqual(sorted(mesh.getGroupsOfNodes()), sorted(mesh.getGroupsOfNodes(True)))
test.assertSequenceEqual(sorted(mesh.getGroupsOfNodes()), sorted(mesh.getGroupsOfNodes(False)))
# do the same thing (compatibily with ParallelMesh)
test.assertTrue(mesh.hasGroupOfNodes("A"))
test.assertTrue(mesh.hasGroupOfNodes("A", True))
test.assertTrue(mesh.hasGroupOfNodes("A", False))
test.assertFalse(mesh.hasGroupOfNodes("AA"))
test.assertFalse(mesh.hasGroupOfNodes("AA", True))
test.assertFalse(mesh.hasGroupOfNodes("AA", False))


# do the same thing (compatibily with ParallelMesh)
test.assertSequenceEqual(sorted(mesh.getGroupsOfCells()), ["Bas", "Haut"])
test.assertSequenceEqual(sorted(mesh.getGroupsOfCells()), sorted(mesh.getGroupsOfCells(False)))
test.assertSequenceEqual(sorted(mesh.getGroupsOfCells()), sorted(mesh.getGroupsOfCells(True)))

# do the same thing (compatibily with ParallelMesh)
test.assertTrue(mesh.hasGroupOfCells("Bas"))
test.assertTrue(mesh.hasGroupOfCells("Bas", True))
test.assertTrue(mesh.hasGroupOfCells("Bas", False))
test.assertFalse(mesh.hasGroupOfCells("Droit"))
test.assertFalse(mesh.hasGroupOfCells("Droit", True))
test.assertFalse(mesh.hasGroupOfCells("Droit", False))

# test coordiantes
coord = mesh.getCoordinates()
test.assertEqual(coord[3], 0.0)
values = coord.getValues()
test.assertEqual(len(values), 27 * 3)

connect = mesh.getConnectivity()
cellsHaut = mesh.getCells("Haut")
test.assertSequenceEqual(cellsHaut, [44, 45, 46, 47])
nodesHaut = mesh.getNodes("Haut")
test.assertSequenceEqual(nodesHaut, [0, 2, 4, 6, 9, 13, 17, 19, 25])
# do the same thing (compatibily with ParallelMesh)
test.assertSequenceEqual(mesh.getNodes("Haut", True), [0, 2, 4, 6, 9, 13, 17, 19, 25])
test.assertSequenceEqual(mesh.getNodes("Haut", True, True), [0, 2, 4, 6, 9, 13, 17, 19, 25])
test.assertSequenceEqual(mesh.getNodes("Haut", False), [0, 2, 4, 6, 9, 13, 17, 19, 25])
test.assertSequenceEqual(mesh.getNodes("Haut", False, False), [0, 2, 4, 6, 9, 13, 17, 19, 25])
test.assertSequenceEqual(mesh.getNodes("Haut", True, False), [0, 2, 4, 6, 9, 13, 17, 19, 25])
test.assertSequenceEqual(mesh.getNodes("Haut", False, True), [0, 2, 4, 6, 9, 13, 17, 19, 25])

# test different variant
test.assertEqual(mesh.getNumberOfNodes(), len(mesh.getNodes()))
test.assertEqual(mesh.getNumberOfCells(), len(mesh.getCells()))

test.assertSequenceEqual(mesh.getNodes(), range(0, mesh.getNumberOfNodes()))
test.assertSequenceEqual(mesh.getCells(), range(0, mesh.getNumberOfCells()))

# do the same thing (compatibily with ParallelMesh)
test.assertSequenceEqual(sorted(mesh.getNodes()), sorted(mesh.getNodes(localNumbering=True)))
test.assertSequenceEqual(
    sorted(mesh.getNodes()), sorted(mesh.getNodes(localNumbering=True, same_rank=True))
)
test.assertSequenceEqual(sorted(mesh.getNodes()), sorted(mesh.getNodes(localNumbering=False)))
test.assertSequenceEqual(
    sorted(mesh.getNodes()), sorted(mesh.getNodes(localNumbering=False, same_rank=True))
)
test.assertSequenceEqual(
    sorted(mesh.getNodes()), sorted(mesh.getNodes(localNumbering=False, same_rank=False))
)
test.assertSequenceEqual(
    sorted(mesh.getNodes()), sorted(mesh.getNodes(localNumbering=True, same_rank=False))
)

medconn = mesh.getMedConnectivity()
medtypes = np.array(mesh.getMedCellsTypes())
test.assertEqual(len(medtypes), mesh.getNumberOfCells())
# cells 1-24: SEG2
test.assertTrue((medtypes[:24] == 102).all())
# cells 25-48: QUAD4
test.assertTrue((medtypes[25:48] == 204).all())
# cells 49-56: HEXA8
test.assertTrue((medtypes[49:] == 308).all())

# check cell #47 (index 46)
quad47 = connect[47 - 1]
test.assertSequenceEqual(quad47, [10, 1, 18, 26])
test.assertSequenceEqual(medconn[47 - 1], [10, 1, 18, 26])
test.assertEqual("M47", mesh.getCellName(47 - 1))
test.assertEqual("N10", mesh.getNodeName(10 - 1))

# always 3 coordinates, even if 'getDimension() == 2'
npcoord = np.array(values).reshape((-1, 3))
# z(cell #47) == 1.
for i in quad47:
    test.assertEqual(npcoord[i - 1][2], 1.0)

# same connectivities for SEG2, QUAD4
test.assertSequenceEqual(connect[12], medconn[12])
test.assertSequenceEqual(connect[36], medconn[36])
# different connectivities for HEXA8
conn_ast = [1, 9, 21, 10, 18, 23, 27, 26]
conn_med = [1, 10, 21, 9, 18, 26, 27, 23]
test.assertSequenceEqual(connect[49 - 1], conn_ast)
test.assertSequenceEqual(medconn[49 - 1], conn_med)

# boundingbox of mesh: (0, 0, 0) -> (1, 1, 1)
test.assertEqual(npcoord.min(), 0.0)
test.assertEqual(npcoord.max(), 1.0)

# test ass new group
mesh.setGroupOfCells("TEST_GMA", [0, 1])
# print(mesh.getGroupsOfCells(), flush=True)
test.assertTrue(mesh.hasGroupOfCells("Haut"))
test.assertTrue(mesh.hasGroupOfCells("TEST_GMA"))
test.assertSequenceEqual(mesh.getCells("TEST_GMA"), [0, 1])
mesh.setGroupOfNodes("TEST_GNO", [14, 8])
test.assertSequenceEqual(mesh.getNodes("TEST_GNO"), [8, 14])

# refine the mesh
mesh = mesh.refine(2)
test.assertEqual(mesh.getNumberOfNodes(), 729)

# read a HEXA27 from ASTER format
mail = code_aster.Mesh()
mail.readAsterFile("zzzz366a.mail")
test.assertTrue(mail.isQuadratic())

m1, m2, m3 = mail.getConnectivity()
# reference connectivities from '.mail' file
ast27 = [
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
]
ast20 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
ast08 = [1, 2, 3, 4, 5, 6, 7, 8]
test.assertSequenceEqual(m1, ast27)
test.assertSequenceEqual(m2, ast20)
test.assertSequenceEqual(m3, ast08)

m1, m2, m3 = mail.getMedConnectivity()
# reference connectivities from IMPR_RESU/MED + mdump
med27 = [
    1,
    4,
    3,
    2,
    5,
    8,
    7,
    6,
    12,
    11,
    10,
    9,
    20,
    19,
    18,
    17,
    13,
    16,
    15,
    14,
    21,
    25,
    24,
    23,
    22,
    26,
    27,
]
med20 = [1, 4, 3, 2, 5, 8, 7, 6, 12, 11, 10, 9, 20, 19, 18, 17, 13, 16, 15, 14]
med08 = [1, 4, 3, 2, 5, 8, 7, 6]
test.assertSequenceEqual(m1, med27)
test.assertSequenceEqual(m2, med20)
test.assertSequenceEqual(m3, med08)

# from ASTER format
mail = code_aster.Mesh()
mail.readAsterFile("ssnp14c.mail")

test.assertFalse(mail.isParallel())
test.assertEqual(mail.getDimension(), 3)
test.assertEqual(mail.getNumberOfNodes(), 22)
test.assertEqual(mail.getNumberOfCells(), 6)
test.assertSequenceEqual(
    mail.getGroupsOfNodes(),
    [
        "NO5",
        "NO6",
        "NO11",
        "NO3",
        "NO16",
        "NO18",
        "NO8",
        "NO12",
        "NO2",
        "NO17",
        "NO20",
        "NO4",
        "NO14",
        "NO7",
        "NO21",
        "NO19",
        "NO15",
        "NO9",
        "NO10",
        "NO1",
    ],
)
test.assertSequenceEqual(
    sorted(mail.getGroupsOfCells()), ["BAS", "DROITE", "GAUCHE", "HAUT", "PENT2"]
)
coord = mail.getCoordinates()
test.assertEqual(coord[3], 1.0)
values = coord.getValues()
test.assertEqual(len(values), 22 * 3)

# refine the mesh
mail = mail.refine(2)
test.assertEqual(mail.getNumberOfNodes(), 505)

# from GMSH format
gmsh = code_aster.Mesh()
gmsh.readGmshFile("ssnv187a.msh")

test.assertFalse(gmsh.isParallel())
test.assertEqual(gmsh.getDimension(), 2)
test.assertEqual(gmsh.getNumberOfNodes(), 132)
test.assertEqual(gmsh.getNumberOfCells(), 207)
test.assertSequenceEqual(gmsh.getGroupsOfNodes(), ["GM5"])
test.assertSequenceEqual(sorted(gmsh.getGroupsOfCells()), ["GM1", "GM3", "GM4", "GM5", "GM6"])
coord = gmsh.getCoordinates()
test.assertEqual(coord[3], 1.0)
values = coord.getValues()
test.assertEqual(len(values), 132 * 3)

# from GIBI format
gibi = code_aster.Mesh()
gibi.readGibiFile("erreu03a.mgib")

test.assertFalse(gibi.isParallel())
test.assertEqual(gibi.getDimension(), 3)
test.assertEqual(gibi.getNumberOfNodes(), 125)
test.assertEqual(gibi.getNumberOfCells(), 84)
test.assertSequenceEqual(sorted(gibi.getGroupsOfNodes()), ["A", "B", "C"])
test.assertSequenceEqual(sorted(gibi.getGroupsOfCells()), ["AB", "BASE1", "CUBE1"])
coord = gibi.getCoordinates()
test.assertEqual(coord[3], 10.0)
values = coord.getValues()
test.assertEqual(len(values), 125 * 3)

# from mesh builder - Cube
builder = code_aster.Mesh.buildCube(refine=3)
test.assertFalse(builder.isParallel())
test.assertEqual(builder.getDimension(), 3)
test.assertEqual(builder.getNumberOfNodes(), 729)
test.assertEqual(builder.getNumberOfCells(), 992)
test.assertSequenceEqual(
    sorted(builder.getGroupsOfNodes()), ["N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8"]
)
test.assertSequenceEqual(
    sorted(builder.getGroupsOfCells()),
    [
        "BACK",
        "BOTTOM",
        "FRONT",
        "LEFT",
        "RIGHT",
        "S13",
        "S21",
        "S24",
        "S26",
        "S34",
        "S37",
        "S51",
        "S56",
        "S68",
        "S75",
        "S78",
        "S84",
        "TOP",
        "VOLUME",
    ],
)
test.assertEqual(code_aster.Mesh.buildCube(refine=2).getNumberOfNodes(), 125)
builder = code_aster.Mesh.buildCube()
test.assertSequenceEqual(builder.getNodesFromCells("VOLUME"), [0, 1, 2, 3, 4, 5, 6, 7])

# from mesh builder - Cylinder
builder = code_aster.Mesh.buildCylinder(refine=3)
test.assertFalse(builder.isParallel())
test.assertEqual(builder.getDimension(), 3)
test.assertEqual(builder.getNumberOfNodes(), 9225)
test.assertEqual(builder.getNumberOfCells(), 10752)
test.assertSequenceEqual(sorted(builder.getGroupsOfNodes()), [])
test.assertSequenceEqual(sorted(builder.getGroupsOfCells()), ["BOTTOM", "SURFEXT", "TOP", "VOLUME"])
test.assertEqual(code_aster.Mesh.buildCylinder(refine=2).getNumberOfNodes(), 1285)
builder = code_aster.Mesh.buildCylinder()
test.assertSequenceEqual(
    builder.getNodesFromCells("BOTTOM"), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
)

# from mesh builder -Square
builder = code_aster.Mesh.buildSquare(refine=3)
test.assertFalse(builder.isParallel())
test.assertEqual(builder.getDimension(), 2)
test.assertEqual(builder.getNumberOfNodes(), 81)
test.assertEqual(builder.getNumberOfCells(), 96)
test.assertSequenceEqual(sorted(builder.getGroupsOfNodes()), ["N1", "N2", "N3", "N4"])
test.assertSequenceEqual(
    sorted(builder.getGroupsOfCells()), ["BOTTOM", "LEFT", "RIGHT", "SURFACE", "TOP"]
)
test.assertSequenceEqual(builder.getNodesFromCells("BOTTOM"), [0, 1, 2, 3, 4, 5, 6, 7, 8])
test.assertEqual(code_aster.Mesh.buildSquare(refine=2).getNumberOfNodes(), 25)

# from mesh builder -  Disk
builder = code_aster.Mesh.buildDisk(refine=3)
test.assertFalse(builder.isParallel())
test.assertEqual(builder.getDimension(), 2)
test.assertEqual(builder.getNumberOfNodes(), 1025)
test.assertEqual(builder.getNumberOfCells(), 1088)
test.assertSequenceEqual(sorted(builder.getGroupsOfNodes()), [])
test.assertSequenceEqual(sorted(builder.getGroupsOfCells()), ["REXT", "SURFACE"])
test.assertSequenceEqual(
    builder.getNodesFromCells("REXT"),
    [
        16,
        32,
        48,
        64,
        80,
        96,
        112,
        128,
        144,
        160,
        176,
        192,
        208,
        224,
        240,
        256,
        272,
        288,
        304,
        320,
        336,
        352,
        368,
        384,
        400,
        416,
        432,
        448,
        464,
        480,
        496,
        512,
        528,
        544,
        560,
        576,
        592,
        608,
        624,
        640,
        656,
        672,
        688,
        704,
        720,
        736,
        752,
        768,
        784,
        800,
        816,
        832,
        848,
        864,
        880,
        896,
        912,
        928,
        944,
        960,
        976,
        992,
        1008,
        1024,
    ],
)
test.assertEqual(code_aster.Mesh.buildDisk(refine=2).getNumberOfNodes(), 257)

# from mesh builder -  Disk
builder = code_aster.Mesh.buildRing(refine=3)
test.assertSequenceEqual(sorted(builder.getGroupsOfCells()), ["REXT", "RINT", "SURFACE"])
test.assertEqual(len(builder.getCells("RINT")), 64)
test.assertEqual(len(builder.getCells("REXT")), 64)

builder = code_aster.Mesh.buildTube(refine=2)
test.assertSequenceEqual(
    sorted(builder.getGroupsOfCells()), ["BOTTOM", "SURFEXT", "SURFINT", "TOP", "VOLUME"]
)
test.assertEqual(len(builder.getCells("SURFEXT")), 128)
test.assertEqual(len(builder.getCells("SURFINT")), 128)

test.printSummary()

code_aster.close()
