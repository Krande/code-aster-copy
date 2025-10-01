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

CA.init("--test", ERREUR=_F(ALARME="ALARME", ERREUR_F="EXCEPTION"))

test = CA.TestCase()

fmt_raison = (
    "-" * 80
    + """

   Exception interceptee
   Message : %s

"""
    + "-" * 80
    + "\n"
)

is_ok = 0
try:
    mesh_3d = LIRE_MAILLAGE(FORMAT="ASTER", UNITE=20)
except CA.AsterError as err:
    print(fmt_raison % str(err))
    # on verifie que l'erreur fatale est bien celle que l'on attendait :
    if err.id_message == "MODELISA4_10":
        is_ok = 1

TAB1 = CREA_TABLE(
    LISTE=(_F(PARA="TEST", TYPE_K="K8", LISTE_K="VALEUR  "), _F(PARA="BOOLEEN", LISTE_I=is_ok))
)
TEST_TABLE(
    REFERENCE="ANALYTIQUE",
    VALE_CALC_I=1,
    VALE_REFE_I=1,
    NOM_PARA="BOOLEEN",
    TABLE=TAB1,
    FILTRE=_F(NOM_PARA="TEST", VALE_K="VALEUR  "),
)

mesh_3d = LIRE_MAILLAGE(FORMAT="ASTER", UNITE=20, VERI_MAIL=_F(VERIF="NON"), INFO=2)
connect = mesh_3d.getConnectivity()

test.assertEqual(mesh_3d.getNumberOfNodes(), 34)
test.assertEqual(mesh_3d.getNumberOfCells(), 29)
test.assertEqual(connect[16], [33, 18, 32])
test.assertEqual(
    connect[28], [9, 8, 10, 11, 4, 5, 7, 6, 28, 32, 29, 31, 23, 21, 25, 27, 16, 19, 18, 17]
)

mesh_3d_fix = mesh_3d.fix(True, info=2)
connect = mesh_3d_fix.getConnectivity()

test.assertEqual(mesh_3d_fix.getNumberOfNodes(), mesh_3d.getNumberOfNodes() - 2)
test.assertEqual(mesh_3d_fix.getNumberOfCells(), mesh_3d.getNumberOfCells())
test.assertEqual(connect[16], [0, 18, 30])
test.assertEqual(
    connect[28], [9, 8, 10, 11, 4, 5, 7, 6, 28, 30, 29, 31, 23, 21, 25, 27, 16, 19, 18, 17]
)

octree_3d = mesh_3d_fix.getOctreeMesh()
test.assertEqual(octree_3d.getNumberOfNodes(), 201)
test.assertEqual(octree_3d.getNumberOfCells(), 92)

is_ok = 0
try:
    mesh_2d = LIRE_MAILLAGE(FORMAT="ASTER", UNITE=21)
except CA.AsterError as err:
    print(fmt_raison % str(err))
    # on verifie que l'erreur fatale est bien celle que l'on attendait :
    if err.id_message == "MODELISA4_10":
        is_ok = 1

TAB1 = CREA_TABLE(
    LISTE=(_F(PARA="TEST", TYPE_K="K8", LISTE_K="VALEUR  "), _F(PARA="BOOLEEN", LISTE_I=is_ok))
)
TEST_TABLE(
    REFERENCE="ANALYTIQUE",
    VALE_CALC_I=1,
    VALE_REFE_I=1,
    NOM_PARA="BOOLEEN",
    TABLE=TAB1,
    FILTRE=_F(NOM_PARA="TEST", VALE_K="VALEUR  "),
)


mesh_2d = LIRE_MAILLAGE(FORMAT="ASTER", UNITE=21, VERI_MAIL=_F(VERIF="NON"), INFO=2)
connect = mesh_2d.getConnectivity()

test.assertEqual(mesh_2d.getNumberOfNodes(), 15)
test.assertEqual(mesh_2d.getNumberOfCells(), 11)
test.assertEqual(connect[8], [14, 1])
test.assertEqual(connect[7], [5, 4, 2, 3, 13, 8, 9, 11])
test.assertEqual(mesh_2d.getNodes("TEST_NO"), [13, 14])
test.assertEqual(mesh_2d.getCells("TEST_MA"), [10])

mesh_2d_fix = mesh_2d.fix(info=2)
connect = mesh_2d_fix.getConnectivity()

test.assertEqual(mesh_2d_fix.getNumberOfNodes(), mesh_2d.getNumberOfNodes() - 2)
test.assertEqual(mesh_2d_fix.getNumberOfCells(), mesh_2d.getNumberOfCells() - 2)
test.assertEqual(connect[8], [0, 1])
test.assertEqual(connect[7], [5, 4, 2, 3, 12, 8, 9, 11])
test.assertEqual(mesh_2d_fix.getNodes("TEST_NO"), [12, 0])
test.assertEqual(mesh_2d_fix.getCells("TEST_MA"), [8])

octree_2d = mesh_2d_fix.getOctreeMesh()
test.assertEqual(octree_2d.getNumberOfNodes(), 35)
test.assertEqual(octree_2d.getNumberOfCells(), 22)

# test for reorientation

mesh_3d = LIRE_MAILLAGE(FORMAT="MED", UNITE=22, INFO=2)
connect = mesh_3d.getConnectivity()

test.assertEqual(connect[0], [5, 1, 8, 16, 20, 17])
test.assertEqual(connect[1], [1, 2, 6, 5, 9, 11, 12, 16])
test.assertEqual(connect[2], [5, 1, 0, 4, 8, 16, 15, 21, 14, 17, 20, 19, 18])

mesh_3d_fix = mesh_3d.fix(True, info=2)

connect = mesh_3d_fix.getConnectivity()

test.assertEqual(connect[0], [5, 8, 1, 17, 20, 16])
test.assertEqual(connect[1], [1, 2, 6, 5, 9, 11, 12, 16])
test.assertEqual(connect[2], [5, 4, 0, 1, 8, 14, 21, 15, 16, 17, 18, 19, 20])

# issue35069
mesh_3d = LIRE_MAILLAGE(FORMAT="MED", UNITE=23, INFO=1)
mesh_3d_fix = mesh_3d.fix(True, info=1)

test.assertEqual(mesh_3d.getNumberOfNodes(), mesh_3d_fix.getNumberOfNodes())
test.assertEqual(mesh_3d.getNumberOfCells(), mesh_3d_fix.getNumberOfCells())

test.printSummary()

FIN()
