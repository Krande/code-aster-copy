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
from code_aster.CA import MPI

test = CA.TestCase()

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

nProc = MPI.ASTER_COMM_WORLD.Get_size()

rank = MPI.ASTER_COMM_WORLD.Get_rank()

pMesh2 = CA.ParallelMesh()
pMesh2.readMedFile("mesh004c/%d.med" % rank, partitioned=True)

model = AFFE_MODELE(
    MAILLAGE=pMesh2,
    AFFE=_F(MODELISATION="D_PLAN", PHENOMENE="MECANIQUE", TOUT="OUI"),
    DISTRIBUTION=_F(METHODE="CENTRALISE"),
)

char_cin = AFFE_CHAR_CINE(
    MODELE=model,
    MECA_IMPO=(
        _F(GROUP_NO="N2", DX=0.0, DY=0.0, DZ=0.0),
        _F(GROUP_NO="N4", DX=0.0, DY=0.0, DZ=0.0),
    ),
)

char_meca = AFFE_CHAR_MECA(
    MODELE=model,
    LIAISON_DDL=_F(GROUP_NO=("N1", "N3"), DDL=("DX", "DX"), COEF_MULT=(1.0, -1.0), COEF_IMPO=0),
    DDL_IMPO=_F(GROUP_NO="N1", DX=1.0),
)

MATER1 = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3))

AFFMAT = AFFE_MATERIAU(MAILLAGE=pMesh2, AFFE=_F(TOUT="OUI", MATER=MATER1))

LI = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=1))

resu = STAT_NON_LINE(
    CHAM_MATER=AFFMAT,
    METHODE="NEWTON",
    COMPORTEMENT=_F(RELATION="ELAS"),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-10),
    EXCIT=(_F(CHARGE=char_cin), _F(CHARGE=char_meca)),
    INCREMENT=_F(LIST_INST=LI),
    MODELE=model,
    NEWTON=_F(MATRICE="TANGENTE", REAC_ITER=1),
    SOLVEUR=_F(METHODE="PETSC", RESI_RELA=1.0e-8, PRE_COND="LDLT_SP"),
)

MyFieldOnNodes = resu.getField("DEPL", 1)
sfon = MyFieldOnNodes.toSimpleFieldOnNodes()

if rank == 0:
    test.assertAlmostEqual(sfon[0, 0], 1.0)
    test.assertAlmostEqual(sfon[2, 0], 0.5175556151165849)
elif rank == 1:
    test.assertAlmostEqual(sfon[0, 0], 1.0)
    test.assertAlmostEqual(sfon[2, 0], 0.5175556151165849)

FIN()
