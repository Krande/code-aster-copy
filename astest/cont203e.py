# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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

from libaster import ContactPairing, ContactComputation, PairingMethod
from code_aster.Commands import *
from code_aster import CA
from code_aster.visu.Pairing import asterPairingProcess as APP
from code_aster.visu.Pairing import pairingObjects as pObj

try:
    import matplotlib

    HAS_MATPLOTLIB_TC = True
except ImportError:
    HAS_MATPLOTLIB_TC = False

DEBUT(CODE="OUI", ERREUR=_F(ALARME="ALARME"))

test = CA.TestCase()
## -------------------------------------------------------
#           Parameters
## -------------------------------------------------------
# - Name of the GROUP_MA for the interface
grma_slv = "CONT_HAUT"
grma_mas = "CONT_BAS"
# - Name of the GROUP_MA for the domain
grma_slv_d = "HAUT"
grma_mas_d = "BAS"
# - Dimension of the mesh
dimMesh = 2


## - Visualisation tests
# - If True, then choose a specific option, else loop over all the options
specific_visu_option = False
# - Option test to run
option_test = 12
# - Dictionnary of options
option_dict = {
    1: {  # - Visualisation of the two solids of the mesh
        "optionMesh": "domain",
        "suboptionMesh": "all",
        "optionPair": "meshOnly",
        "addMeshNodes": True,
        "addLegend": True,
        "index": None,
        "indexPlaneProjected": None,
    },
    2: {  # - Visualisation of the two interfaces of the mesh
        "optionMesh": "interface",
        "suboptionMesh": "all",
        "optionPair": "meshOnly",
        "addMeshNodes": True,
        "addLegend": True,
        "index": None,
        "indexPlaneProjected": None,
    },
    3: {  # - Visualisation of all the pairs on the two interfaces of the mesh
        "optionMesh": "interface",
        "suboptionMesh": "all",
        "optionPair": "pairs",
        "addMeshNodes": True,
        "addLegend": True,
        "index": None,
        "indexPlaneProjected": None,
    },
    4: {  # - Visualisation of the pair of index = 0 (first pair) on the two interfaces of the mesh
        "optionMesh": "interface",
        "suboptionMesh": "givenPair",
        "optionPair": "pairs",
        "addMeshNodes": True,
        "addLegend": True,
        "index": 0,
        "indexPlaneProjected": None,
    },
    5: {  # - Visualisation of the pair which involves the slave cell of index = 17
        #     on the two interfaces of the mesh
        "optionMesh": "interface",
        "suboptionMesh": "givenSlvIndex",
        "optionPair": "pairs",
        "addMeshNodes": True,
        "addLegend": True,
        "index": 17,
        "indexPlaneProjected": None,
    },
    6: {  # - Visualisation of the intersection points
        # of index = 0 (first pair) on the two interfaces of the mesh
        "optionMesh": "interface",
        "suboptionMesh": "givenPair",
        "optionPair": "intePoints",
        "addMeshNodes": True,
        "addLegend": True,
        "index": 0,
        "indexPlaneProjected": None,
    },
    7: {  # - Visualisation of the intersection points of the pairs
        #     which involve the slave cell of index = 17
        #     on the two interfaces of the mesh
        "optionMesh": "interface",
        "suboptionMesh": "givenSlvIndex",
        "optionPair": "intePoints",
        "addMeshNodes": True,
        "addLegend": True,
        "index": 17,
        "indexPlaneProjected": None,
    },
    8: {  # - Visualisation of the quadrature points
        # of index = 0 (first pair) on the two interfaces of the mesh
        "optionMesh": "interface",
        "suboptionMesh": "givenPair",
        "optionPair": "quadPoints",
        "addMeshNodes": True,
        "addLegend": True,
        "index": 0,
        "indexPlaneProjected": None,
    },
    9: {  # - Visualisation of the quadrature points of the pairs
        #     which involve the slave cell of index = 17
        #     on the two interfaces of the mesh
        "optionMesh": "interface",
        "suboptionMesh": "givenSlvIndex",
        "optionPair": "quadPoints",
        "addMeshNodes": True,
        "addLegend": True,
        "index": 17,
        "indexPlaneProjected": None,
    },
    10: {  # - Visualisation of all the cells paired to the slave cell of index = 18
        "optionMesh": "selectSlvCell",
        "suboptionMesh": "givenSlvIndex",
        "optionPair": "pairs",
        "addMeshNodes": True,
        "addLegend": True,
        "index": 18,
        "indexPlaneProjected": None,
    },
    11: {  # - Visualisation of all the intersections of
        # the pairs to which the slave cell of index = 18 belongs
        "optionMesh": "selectSlvCell",
        "suboptionMesh": "givenSlvIndex",
        "optionPair": "intePoints",
        "addMeshNodes": True,
        "addLegend": True,
        "index": 18,
        "indexPlaneProjected": None,
    },
    12: {  # - Visualisation of all the quadrature points of
        # the pairs to which the slave cell of index = 18 belongs
        "optionMesh": "selectSlvCell",
        "suboptionMesh": "givenSlvIndex",
        "optionPair": "quadPoints",
        "addMeshNodes": True,
        "addLegend": True,
        "index": 18,
        "indexPlaneProjected": None,
    },
}

## -------------------------------------------------------
#           Definition of the mesh and the model
## -------------------------------------------------------
ma = LIRE_MAILLAGE(FORMAT="MED")

ma = MODI_MAILLAGE(
    reuse=ma, MAILLAGE=ma, ORIE_PEAU_3D=(_F(GROUP_MA=grma_slv), _F(GROUP_MA=grma_mas))
)

modele = AFFE_MODELE(
    MAILLAGE=ma, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION=("D_PLAN",))
)

## -------------------------------------------------------
#           Get aster information
## -------------------------------------------------------
# - Check initialisation procedure
AsterPairing = APP.AsterPairingProcess(grma_slv, grma_mas, ma)
AsterPairing.setMethod(method="BrutForce")

test.assertEqual(AsterPairing._groupMaSlv, grma_slv, msg="Init. Slave group name")
test.assertEqual(AsterPairing._groupMaMas, grma_mas, msg="Init. Master group name")
test.assertEqual(AsterPairing._method, "BrutForce", msg="Init. Pairing algo")
test.assertEqual(AsterPairing._asterMesh, ma, msg="Init. aster mesh")

# - Pairing procedure
AsterPairing.run()
test.assertEqual(len(AsterPairing._listPairs), 7, msg="Check pairing - number of pairs")
test.assertEqual(
    len(AsterPairing._listPairs),
    len(AsterPairing._intePointsList),
    msg="Check pairing - consistency sizes of list (1)",
)
test.assertEqual(
    len(AsterPairing._listPairs),
    len(AsterPairing._quadPointsList),
    msg="Check pairing - consistency sizes of list (2)",
)

# - Connectivity and group information (for visualisation)
AsterPairing.extractMeshInfosFromAsterMesh()
test.assertEqual(AsterPairing._coords.shape[0], 52, msg="Check connectivity array (1)")
test.assertEqual(AsterPairing._coords.shape[1], 3, msg="Check connectivity array (2)")
test.assertEqual(len(AsterPairing._asterConnectivity), 66, msg="Check connectivity array (3)")

## -------------------------------------------------------
#           Prepare visualisation
## -------------------------------------------------------
# - PairingAnalysisAster class
AsterVisu = pObj.PairingAnalysisAster(
    dimMesh, grma_mas_d, grma_mas, grma_slv_d, grma_slv, AsterPairing
)

test.assertTrue(AsterVisu._flag_MeshInfos, msg="Check PairingAnalysisAster flag (1)")
test.assertTrue(AsterVisu._flag_PairingInfos, msg="Check PairingAnalysisAster flag (2)")
test.assertIsNotNone(AsterVisu._listPairs, msg="Check PairingAnalysisAster flag (3)")
test.assertIsNotNone(AsterVisu._listIntersectionPts, msg="Check PairingAnalysisAster flag (4)")
test.assertIsNotNone(AsterVisu._listQuadraturePts, msg="Check PairingAnalysisAster flag (5)")

## -------------------------------------------------------
#           Plot interactive
## -------------------------------------------------------
if HAS_MATPLOTLIB_TC:
    if specific_visu_option:
        AsterVisu.plotInteractive(
            option_dict[option_test]["optionMesh"],
            option_dict[option_test]["suboptionMesh"],
            option_dict[option_test]["optionPair"],
            addMeshNodes=option_dict[option_test]["addMeshNodes"],
            addLegend=option_dict[option_test]["addLegend"],
            index=option_dict[option_test]["index"],
            indexPlaneProjected=option_dict[option_test]["indexPlaneProjected"],
        )
    else:
        for option_visu in list(option_dict.keys()):
            AsterVisu.plotInteractive(
                option_dict[option_visu]["optionMesh"],
                option_dict[option_visu]["suboptionMesh"],
                option_dict[option_visu]["optionPair"],
                addMeshNodes=option_dict[option_visu]["addMeshNodes"],
                addLegend=option_dict[option_visu]["addLegend"],
                index=option_dict[option_visu]["index"],
                indexPlaneProjected=option_dict[option_visu]["indexPlaneProjected"],
            )
FIN()
