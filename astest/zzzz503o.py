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
# Test to validate IMPR_RESU FICHIER_UNIQUE with NOM_CMP and new partitioning

import os.path as osp


from code_aster.CA import MPI
from code_aster.Commands import *
from code_aster import CA
from code_aster.Utilities import SharedTmpdir

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

rank = MPI.ASTER_COMM_WORLD.Get_rank()

test = CA.TestCase()

# --------------------------------------------------------------------------------------------------
# Valeur de chargement de reference
# --------------------------------------------------------------------------------------------------
q_mdr = 325
q_elv = 380
# --------------------------------------------------------------------------------------------------


# MA = LIRE_MAILLAGE()


def checkJoints(mesh):
    comm = MPI.ASTER_COMM_WORLD
    l2G = mesh.getLocalToGlobalNodeIds()

    j = 0
    for proc in mesh.getOppositeDomains():
        fJ = mesh.getSendJoint(j)
        gFJ = []
        for i in range(int(len(fJ) / 2)):
            gFJ.append(l2G[fJ[2 * i] - 1])

        sJ = mesh.getReceiveJoint(j)
        gSJ = []
        for i in range(int(len(sJ) / 2)):
            gSJ.append(l2G[sJ[2 * i] - 1])

        if proc < rank:
            comm.send(gFJ, dest=proc, tag=j)
            data1 = comm.recv(source=proc, tag=j)
            if not data1 == gFJ:
                break
            test.assertEqual(data1 == gFJ, True)
        else:
            data1 = comm.recv(source=proc, tag=j)
            comm.send(gSJ, dest=proc, tag=j)
            if not data1 == gSJ:
                break
            test.assertEqual(data1 == gSJ, True)
        j += 1


mesh3 = CA.IncompleteMesh()
mesh3.readMedFile("zzzz503o.mmed")
bMesh = CA.MeshBalancer()
bMesh.buildFromBaseMesh(mesh3)
meshGraph = CA.MeshConnectionGraph()
meshGraph.buildFromIncompleteMesh(mesh3)
part2 = CA.PtScotchPartitioner()
part2.buildGraph(meshGraph)
scotchPart = part2.partitionGraph()
MA = bMesh.applyBalancingStrategy(scotchPart)
checkJoints(MA)

# MA = CA.ParallelMesh()

MA = CREA_MAILLAGE(MAILLAGE=MA, QUAD_LINE=_F(TOUT="OUI"))

MO = AFFE_MODELE(
    MAILLAGE=MA,
    AFFE=(
        _F(
            GROUP_MA=("TUYAU", "SECTION_0", "SECTION_1", "INTRADOS"),
            PHENOMENE="MECANIQUE",
            MODELISATION="3D",
        ),
        _F(GROUP_MA=("POINT_0", "POINT_1"), PHENOMENE="MECANIQUE", MODELISATION="DIS_TR"),
    ),
)


# =================================================================================================================
#   CARACTERISTIQUES MECANIQUES
# =================================================================================================================


A42P1 = DEFI_MATERIAU(ELAS=_F(E=204000.0, NU=0.3), ECRO_LINE=_F(SY=215, D_SIGM_EPSI=2000))


CHMAT = AFFE_MATERIAU(MAILLAGE=MA, AFFE=(_F(GROUP_MA="TUYAU", MATER=A42P1),))


CHCARA = AFFE_CARA_ELEM(
    MODELE=MO, DISCRET=_F(GROUP_MA=("POINT_0", "POINT_1"), CARA="K_TR_D_N", VALE=(0, 0, 0, 0, 0, 0))
)


# =================================================================================================================
#   LIAISONS ET CHARGEMENTS
# =================================================================================================================


LIAISON0 = AFFE_CHAR_MECA(
    MODELE=MO,
    LIAISON_ELEM=(
        _F(OPTION="3D_POU", GROUP_MA_1="SECTION_0", GROUP_MA_2="POINT_0"),
        _F(OPTION="3D_POU", GROUP_MA_1="SECTION_1", GROUP_MA_2="POINT_1"),
    ),
)

LIAISON1 = AFFE_CHAR_CINE(
    MODELE=MO, MECA_IMPO=_F(GROUP_MA="POINT_0", DX=0, DY=0, DZ=0, DRX=0, DRY=0, DRZ=0)
)

LIAISON2 = AFFE_CHAR_CINE(MODELE=MO, MECA_IMPO=_F(GROUP_MA="SECTION_1", DZ=0))


EFFORT = AFFE_CHAR_MECA(
    MODELE=MO, FORCE_NODALE=(_F(GROUP_NO="POINT_1", FX=1.36394e1, FY=-2.80316e1, MZ=2.18814e5),)
)


PRESSION = AFFE_CHAR_MECA(MODELE=MO, PRES_REP=_F(GROUP_MA="INTRADOS", PRES=2.35009e-2))


# =================================================================================================================
#   DEROULEMENT DU CALCUL
# =================================================================================================================

DISCRET = DEFI_LIST_INST(
    DEFI_LIST=_F(VALE=(0, 200, 225, 250, 275, 300, q_mdr, 350, 360, 370, q_elv)),
    ECHEC=_F(SUBD_NIVEAU=2, SUBD_PAS=2),
)

DISCRET = DEFI_LIST_REEL(VALE=(0, 200, 225, 250, 275, 300, q_mdr, 350, 360, 370, q_elv))


RAMPE = DEFI_FONCTION(NOM_PARA="INST", VALE=(0, 0, 1, 1), PROL_DROITE="LINEAIRE")


snl = dict(
    MODELE=MO,
    CARA_ELEM=CHCARA,
    CHAM_MATER=CHMAT,
    EXCIT=(
        _F(CHARGE=LIAISON0),
        _F(CHARGE=LIAISON1),
        _F(CHARGE=LIAISON2),
        _F(CHARGE=EFFORT, FONC_MULT=RAMPE),
        _F(CHARGE=PRESSION, FONC_MULT=RAMPE),
    ),
    COMPORTEMENT=(
        _F(TOUT="OUI", RELATION="ELAS", DEFORMATION="PETIT"),
        _F(GROUP_MA="TUYAU", RELATION="ELAS", DEFORMATION="PETIT"),
    ),
    CONVERGENCE=_F(
        ITER_GLOB_MAXI=10,
        RESI_REFE_RELA=1.0e-3,
        SIGM_REFE=215,
        EFFORT_REFE=215 * 2 * 2,
        MOMENT_REFE=215 * 2 * 2 * 100,
    ),
    NEWTON=_F(PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    SOLVEUR=_F(
        METHODE="PETSC",
        PRE_COND="SANS",
        OPTION_PETSC="-pc_type lu -pc_factor_mat_solver_type mumps  -ksp_monitor_true_residual",
    ),
    INFO=1,
)


EVOL = STAT_NON_LINE(INCREMENT=_F(LIST_INST=DISCRET, INST_FIN=200), **snl)

with SharedTmpdir("zzzz503d_") as tmpdir:
    medfile = osp.join(tmpdir.path, "resu.zzzz503d.no.med")
    DEFI_FICHIER(UNITE=99, FICHIER=medfile, TYPE="BINARY")
    IMPR_RESU(
        FORMAT="MED",
        FICHIER_UNIQUE="OUI",
        INFO=2,
        UNITE=99,
        RESU=(
            _F(
                GROUP_MA="TUYAU",
                INST=200,
                RESULTAT=EVOL,
                NOM_CHAM="DEPL",
                NOM_CHAM_MED="DEPL",
                NOM_CMP=("DX", "DY", "DZ"),
            ),
        ),
    )
    DEFI_FICHIER(ACTION="LIBERER", UNITE=99)

IMPR_RESU(
    FORMAT="MED",
    UNITE=61,
    RESU=(
        _F(
            GROUP_MA="TUYAU",
            INST=200,
            RESULTAT=EVOL,
            NOM_CHAM="DEPL",
            NOM_CHAM_MED="DEPL",
            NOM_CMP=("DX", "DY", "DZ"),
        ),
    ),
)

FIN()
