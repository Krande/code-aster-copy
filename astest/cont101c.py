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

# -- Lines for dummy test and visualisation test
path_repo_save = (
    "/home/i48378/Documents/workspace/sandbox/sandbox-contact-mec-pairing-tests/cont101c"
)
import os.path as osp

# case - a (Slave side = CONT_BAS)
path_pairs_a = osp.join(path_repo_save, "case_a/listOfPairs.npy")
path_intepoints_a = osp.join(path_repo_save, "case_a/intersectionPoints.npy")
path_quadpoints_a = osp.join(path_repo_save, "case_a/quadraturePoints.npy")
# case - b (Slave side = CONT_HAUT)
path_pairs_b = osp.join(path_repo_save, "case_b/listOfPairs.npy")
path_intepoints_b = osp.join(path_repo_save, "case_b/intersectionPoints.npy")
path_quadpoints_b = osp.join(path_repo_save, "case_b/quadraturePoints.npy")
save_data = True

# -- Core of the test case
from code_aster.Commands import *
from code_aster import CA
from libaster import PairingMethod


DEBUT(
    CODE="OUI",
    ERREUR=_F(ALARME="ALARME"),
    # DEBUG=_F(SDVERI='OUI',),
    INFO=1,
)

test = CA.TestCase()

Mail = LIRE_MAILLAGE(UNITE=20, FORMAT="MED", INFO=2)

Mail = MODI_MAILLAGE(
    reuse=Mail, MAILLAGE=Mail, ORIE_PEAU=_F(GROUP_MA_PEAU=("CONT_HAUT", "CONT_BAS"))
)

MODI = AFFE_MODELE(MAILLAGE=Mail, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

# # Slave side - CONT_BAS
# DEFICO_BAS = DEFI_CONT(
#     MODELE=MODI,
#     INFO=2,
#     ZONE=(
#         _F(
#             APPARIEMENT="MORTAR",
#             GROUP_MA_MAIT="CONT_HAUT",
#             GROUP_MA_ESCL="CONT_BAS",
#             ALGO_CONT="LAGRANGIEN",
#             CONTACT_INIT="OUI",
#         ),
#     ),
# )

# # Pairing checks
# pair = ContactPairing(DEFICO_BAS)
# pair.compute()

# # some checks for the zone 0:
# zone = DEFICO_BAS.getContactZone(0)
# meshPair = zone.getMeshPairing()
# listPairs_zone = meshPair.getListOfPairs()

# Generate pairs
meshPair = CA.MeshPairing()
meshPair.setMesh(Mail)
meshPair.setVerbosity(2)
meshPair.setPair("CONT_BAS", "CONT_HAUT")
meshPair.setMethod(PairingMethod.Fast)
meshPair.compute()

# Get pairs
nbPairs = meshPair.getNumberOfPairs()
listPairs = meshPair.getListOfPairs()

print("Nb pairs: ", nbPairs)
print("List pairs: ", listPairs)


# Some checks
test.assertSequenceEqual(
    listPairs,
    [
        (186, 319),
        (198, 319),
        (198, 314),
        (198, 309),
        (187, 309),
        (202, 309),
        (202, 314),
        (202, 310),
        (202, 319),
        (202, 315),
        (202, 320),
        (202, 316),
        (189, 309),
        (189, 304),
        (203, 316),
        (203, 315),
        (203, 311),
        (203, 310),
        (203, 306),
        (203, 309),
        (203, 305),
        (197, 316),
        (197, 315),
        (197, 321),
        (197, 320),
        (197, 319),
        (188, 304),
        (201, 304),
        (201, 309),
        (201, 305),
        (201, 310),
        (201, 306),
        (195, 306),
        (195, 311),
        (195, 316),
        (195, 312),
        (195, 317),
        (196, 321),
        (196, 322),
        (196, 316),
        (190, 304),
        (200, 306),
        (200, 305),
        (200, 304),
        (193, 317),
        (193, 312),
        (193, 311),
        (193, 307),
        (193, 306),
        (194, 317),
        (194, 316),
        (194, 322),
        (194, 321),
        (191, 304),
        (199, 306),
        (199, 307),
        (192, 307),
    ],
)

intePoints = meshPair.getIntersectionPoints(0)

print("intePoints: ", intePoints)

test.assertAlmostEqual(intePoints[0][0], 50.0)
test.assertAlmostEqual(intePoints[0][1], 183.3333333333333)
test.assertAlmostEqual(intePoints[0][2], 200.0)


# if save_data:
#     np.save(path_pairs_a, listPairs)
#     intePointsList = []
#     quadPointsList = []
#     # - Loop over the pairs
#     for index_pair in range(nbPairs):
#         # - Store and save intersection points
#         IntePoints = meshPair.getIntersectionPoints(index_pair)
#         intePts_current = [tuple(intePt) for intePt in IntePoints]
#         intePointsList.append(intePts_current)
#         # - Store and save integration points
#         quadPoints = meshPair.getQuadraturePoints(index_pair)
#         quadPts_current = [tuple(quaPt) for quaPt in quadPoints]
#         quadPointsList.append(quadPts_current)
#     # - Save data
#     np.save(path_intepoints_a, intePointsList)
#     np.save(path_quadpoints_a, quadPointsList)

# # Slave side - CONT_BAS
# DEFICO_HAUT = DEFI_CONT(
#     MODELE=MODI,
#     INFO=2,
#     ZONE=(
#         _F(
#             APPARIEMENT="MORTAR",
#             GROUP_MA_MAIT="CONT_BAS",
#             GROUP_MA_ESCL="CONT_HAUT",
#             ALGO_CONT="LAGRANGIEN",
#             CONTACT_INIT="OUI",
#         ),
#     ),
# )

# # Pairing checks
# pair = ContactPairing(DEFICO_HAUT)
# pair.compute()

# # some checks for the zone 0:
# zone = DEFICO_HAUT.getContactZone(0)
# meshPair = zone.getMeshPairing()
# listPairs_zone = meshPair.getListOfPairs()

# # Generate pairs
# meshPair = CA.MeshPairing(Mail)
# meshPair.setVerbosity(2)
# meshPair.setPair("CONT_HAUT", "CONT_BAS")
# meshPair.setMethod("NEW")
# meshPair.compute()

# # Get pairs
# nbPairs = meshPair.getNumberOfPairs()
# listPairs = meshPair.getListOfPairs()

# # Some checks
# test.assertSequenceEqual(listPairs,
#         [(304, 188), (304, 190), (304, 189), (304, 191), (304, 201), (304, 200), (309, 201), (309, 189), (309, 203), \
#         (309, 187), (309, 202), (309, 198), (305, 200), (305, 201), (305, 203), (314, 198), (314, 202), (310, 202), \
#         (310, 203), (310, 201), (306, 203), (306, 201), (306, 195), (306, 200), (306, 193), (306, 199), (319, 202), \
#         (319, 198), (319, 197), (319, 186), (315, 202), (315, 203), (315, 197), (311, 203), (311, 195), (311, 193), \
#         (307, 199), (307, 192), (307, 193), (320, 197), (320, 202), (316, 197), (316, 202), (316, 196), (316, 203), \
#         (316, 194), (316, 195), (312, 193), (312, 195), (321, 197), (321, 196), (321, 194), (317, 195), (317, 193), \
#         (317, 194), (322, 194), (322, 196)])

# if save_data:
#     intePointsList = []
#     quadPointsList = []
#     # - Loop over the pairs
#     for index_pair in range(nbPairs):
#         # - Store and save intersection points
#         IntePoints = meshPair.getIntersectionPoints(index_pair)
#         intePts_current = [tuple(intePt) for intePt in IntePoints]
#         intePointsList.append(intePts_current)
#         # - Store and save integration points
#         quadPoints = meshPair.getQuadraturePoints(indexPair)
#         quadPts_current = [tuple(quaPt) for quaPt in quadPoints]
#         quadPointsList.append(quadPts_current)
#     # - Save data
#     np.save(path_intepoints_b, intePointsList)
#     np.save(path_quadpoints_b, quadPointsList)


FIN()
