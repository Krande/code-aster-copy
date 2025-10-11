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

# This a verification test to verify that AFFE_CHAR_MECA EFFE_FOND
# works fine in parallel

import numpy as np
from code_aster.Commands import *
from code_aster import CA
from code_aster.CA import MPI

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"))

test = CA.TestCase()


def buildCompleteFieldOnNodes(field):
    """
    Build complete (over processes) field on nodes
    Arguments:
        field (SimpleFieldOnNodes): field to complete

    Returns:
        list: list of lists containing all field values over processes on global numbering
    """
    field.updateValuePointers()
    mesh = field.getMesh()
    lTGN = mesh.getLocalToGlobalNodeIds()
    maxNodes = max(lTGN)
    maxNodes = MPI.ASTER_COMM_WORLD.allreduce(maxNodes, MPI.MAX)

    innerNodesSet = set(mesh.getInnerNodes())
    values, mask = field.getValues()
    nbNode = field.getNumberOfNodes()
    nbCmp = field.getNumberOfComponents()

    completeField = np.zeros((maxNodes + 1, nbCmp))
    cmpt = 0
    for idNode in range(nbNode):
        if idNode in innerNodesSet:
            globNodeId = lTGN[idNode]
            toAdd = []
            for iCmp in range(nbCmp):
                toAdd.append(values[idNode, iCmp])
            completeField[globNodeId] = np.array(toAdd)
        cmpt += 1

    return MPI.ASTER_COMM_WORLD.allreduce(completeField, MPI.SUM)


# First: sequential reference
M = CA.Mesh()
M.readMedFile("fort.20")

MO = AFFE_MODELE(MAILLAGE=M, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

MODI_MAILLAGE(reuse=M, MAILLAGE=M, ORIE_PEAU=_F(GROUP_MA_PEAU="Bas"))

MA = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3, ALPHA=0.0))

CM = AFFE_MATERIAU(MAILLAGE=M, AFFE=_F(TOUT="OUI", MATER=MA))

CH1 = AFFE_CHAR_CINE(MODELE=MO, INFO=2, MECA_IMPO=(_F(GROUP_NO="ABCD", DX=0.0, DY=0.0, DZ=0.0),))

CH2 = AFFE_CHAR_MECA(
    MODELE=MO,
    INFO=2,
    PRES_REP=_F(GROUP_MA="Interieur", PRES=60000.0),
    EFFE_FOND=(_F(GROUP_MA_INT="Contour", GROUP_MA="Bas", PRES=60000.0),),
)

RESU = MECA_STATIQUE(MODELE=MO, CHAM_MATER=CM, EXCIT=(_F(CHARGE=CH1), _F(CHARGE=CH2)))

depl = RESU.getField("DEPL", 1)
flatDeplRef = depl.getValues()

rank = CA.MPI.ASTER_COMM_WORLD.Get_rank()

# First: parallel comparison
mesh = CA.IncompleteMesh()
mesh.readMedFile("fort.20", verbose=0)
bMesh = CA.MeshBalancer()
bMesh.buildFromBaseMesh(mesh)
M = None
listP0 = [1, 2, 3, 4, 9, 10, 15, 16]
listP1 = [5, 6, 7, 8, 11, 12, 13, 14]
if rank == 0:
    M = bMesh.applyBalancingStrategy(listP0)
elif rank == 1:
    M = bMesh.applyBalancingStrategy(listP1)

MO = AFFE_MODELE(MAILLAGE=M, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

MODI_MAILLAGE(reuse=M, MAILLAGE=M, ORIE_PEAU=_F(GROUP_MA_PEAU="Bas"))

MA = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3, ALPHA=0.0))

CM = AFFE_MATERIAU(MAILLAGE=M, AFFE=_F(TOUT="OUI", MATER=MA))

CH1 = AFFE_CHAR_CINE(MODELE=MO, INFO=2, MECA_IMPO=(_F(GROUP_NO="ABCD", DX=0.0, DY=0.0, DZ=0.0),))

CH2 = AFFE_CHAR_MECA(
    MODELE=MO,
    INFO=2,
    PRES_REP=_F(GROUP_MA="Interieur", PRES=60000.0),
    EFFE_FOND=(_F(GROUP_MA_INT="Contour", GROUP_MA="Bas", PRES=60000.0),),
)

RESU = MECA_STATIQUE(MODELE=MO, CHAM_MATER=CM, EXCIT=(_F(CHARGE=CH1), _F(CHARGE=CH2)))

depl = RESU.getField("DEPL", 1)
compDepl = buildCompleteFieldOnNodes(depl.toSimpleFieldOnNodes())
flatDepl = [x for xs in compDepl for x in xs]

for init, comp in zip(flatDeplRef, flatDepl):
    test.assertAlmostEqual(init, comp, delta=1e-13)

FIN()
