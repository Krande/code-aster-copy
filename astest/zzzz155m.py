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

from code_aster.CA import MPI
from code_aster.Commands import *
from code_aster import CA

# Convergence verification in STAT_NON_LINE with LIAISON_MAIL
# Sequential/Parallel comparison

rank = MPI.ASTER_COMM_WORLD.Get_rank()
size = MPI.ASTER_COMM_WORLD.Get_size()

import numpy as np

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

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


# First: Parallel
filename = "zzzz155m.med"
mesh = CA.ParallelMesh()
mesh.readMedFile(filename)

MAT = DEFI_MATERIAU(ELAS=_F(E=9.84e4, NU=0.3))

MATC = DEFI_MATERIAU(DIS_CONTACT=_F(RIGI_NOR=1e6, RIGI_TAN=1e6, COULOMB=0.0))

MODE = AFFE_MODELE(
    MAILLAGE=mesh, AFFE=(_F(GROUP_MA=("vol1", "vol2"), PHENOMENE="MECANIQUE", MODELISATION="3D"),)
)

MATE = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=(_F(GROUP_MA=("vol1", "vol2"), MATER=MAT),))

BLOC = AFFE_CHAR_CINE(MODELE=MODE, MECA_IMPO=(_F(GROUP_MA=("gauche"), DX=0, DY=0, DZ=0),))

CHAR = AFFE_CHAR_CINE(MODELE=MODE, MECA_IMPO=(_F(GROUP_MA=("droite"), DX=1, DY=1),))

COLL = AFFE_CHAR_MECA(
    MODELE=MODE, LIAISON_MAIL=(_F(GROUP_NO_ESCL="collage", GROUP_MA_MAIT="vol2"),)
)

LINST = DEFI_LIST_REEL(DEBUT=0, INTERVALLE=(_F(JUSQU_A=2, NOMBRE=3),))

LINST2 = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=LINST), ECHEC=_F(SUBD_NIVEAU=5, SUBD_PAS=10))

LINEDEPL = DEFI_FONCTION(NOM_PARA="INST", ABSCISSE=(0, 1, 2), ORDONNEE=(0, 8, 0))

RESU = STAT_NON_LINE(
    MODELE=MODE,
    CHAM_MATER=MATE,
    EXCIT=(_F(CHARGE=BLOC), _F(CHARGE=COLL), _F(CHARGE=CHAR, FONC_MULT=LINEDEPL)),
    INCREMENT=_F(LIST_INST=LINST2),
    COMPORTEMENT=(_F(DEFORMATION="PETIT", RELATION="ELAS", GROUP_MA=("vol1", "vol2")),),
    CONVERGENCE=_F(ITER_GLOB_MAXI=100, RESI_GLOB_RELA=1e-6),
    SOLVEUR=_F(METHODE="PETSC", PRE_COND="LDLT_DP"),
)

# Then: Sequential
mesh2 = CA.Mesh()
mesh2.readMedFile(filename)

MOD2 = AFFE_MODELE(
    MAILLAGE=mesh2,
    DISTRIBUTION=_F(METHODE="CENTRALISE"),
    AFFE=(_F(GROUP_MA=("vol1", "vol2"), PHENOMENE="MECANIQUE", MODELISATION="3D"),),
)

MAT2 = AFFE_MATERIAU(MAILLAGE=mesh2, AFFE=(_F(GROUP_MA=("vol1", "vol2"), MATER=MAT),))

BLO2 = AFFE_CHAR_CINE(MODELE=MOD2, MECA_IMPO=(_F(GROUP_MA=("gauche"), DX=0, DY=0, DZ=0),))

CHA2 = AFFE_CHAR_CINE(MODELE=MOD2, MECA_IMPO=(_F(GROUP_MA=("droite"), DX=1, DY=1),))

COL2 = AFFE_CHAR_MECA(
    MODELE=MOD2, LIAISON_MAIL=(_F(GROUP_NO_ESCL="collage", GROUP_MA_MAIT="vol2"),)
)

RES2 = STAT_NON_LINE(
    MODELE=MOD2,
    CHAM_MATER=MAT2,
    EXCIT=(_F(CHARGE=BLO2), _F(CHARGE=COL2), _F(CHARGE=CHA2, FONC_MULT=LINEDEPL)),
    INCREMENT=_F(LIST_INST=LINST2),
    COMPORTEMENT=(_F(DEFORMATION="PETIT", RELATION="ELAS", GROUP_MA=("vol1", "vol2")),),
    CONVERGENCE=_F(ITER_GLOB_MAXI=100, RESI_GLOB_RELA=1e-6),
    SOLVEUR=_F(METHODE="PETSC", PRE_COND="LDLT_DP"),
)

incompleteField1 = RESU.getField("DEPL", 2).toSimpleFieldOnNodes()
# Parallel field completion
field1 = buildCompleteFieldOnNodes(incompleteField1)

field2 = RES2.getField("DEPL", 2).toSimpleFieldOnNodes()

values1 = field1
values2 = field2.getValues()[0]

# 1e-14 sequential/parallel comparison after 2 iterations
for initArray, newArray in zip(values1, values2):
    for initVal, newVal in zip(initArray, newArray):
        test.assertAlmostEqual(initVal, newVal, 14)

FIN()
