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

rank = MPI.ASTER_COMM_WORLD.Get_rank()

pMesh2 = CA.ParallelMesh()
pMesh2.readMedFile("zzzz504f/%d.med" % rank, partitioned=True)

model = AFFE_MODELE(
    MAILLAGE=pMesh2, AFFE=_F(MODELISATION="D_PLAN", PHENOMENE="MECANIQUE", TOUT="OUI")
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
    LIAISON_RBE3=_F(
        GROUP_NO_MAIT="N5",
        GROUP_NO_ESCL=("N1", "N2", "N3", "N4"),
        COEF_ESCL=1,
        DDL_MAIT="DX",
        DDL_ESCL="DX-DY-DZ",
    ),
    DDL_IMPO=_F(GROUP_NO="N5", DX=0.5),
)
char_meca.debugPrint(10 + rank)

MATER1 = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3))

AFFMAT = AFFE_MATERIAU(MAILLAGE=pMesh2, AFFE=_F(TOUT="OUI", MATER=MATER1))

resu = MECA_STATIQUE(
    CHAM_MATER=AFFMAT,
    MODELE=model,
    EXCIT=(_F(CHARGE=char_cin), _F(CHARGE=char_meca)),
    SOLVEUR=_F(METHODE="PETSC", PRE_COND="SANS", RESI_RELA=1.0e-10),
)
resu.debugPrint(10 + rank)

resu.printMedFile("test" + str(rank) + ".med")
# from shutil import copyfile
# copyfile("test"+str(rank)+".med", "/home/siavelis/test"+str(rank)+".med")

MyFieldOnNodes = resu.getField("DEPL", 1)
sfon = MyFieldOnNodes.toSimpleFieldOnNodes()
sfon.debugPrint(10 + rank)
sfon.build()

# DX displacement on nodes "N1" and "N3", comparison with sequential results
TEST_RESU(
    RESU=(
        _F(
            GROUP_NO=("N1",),
            NOM_CHAM="DEPL",
            NOM_CMP="DX",
            NUME_ORDRE=1,
            RESULTAT=resu,
            VALE_CALC=0.979591885789513,
        ),
        _F(
            GROUP_NO=("N3",),
            NOM_CHAM="DEPL",
            NOM_CMP="DX",
            NUME_ORDRE=1,
            RESULTAT=resu,
            VALE_CALC=1.020408114214617,
        ),
    )
)

RAMPE = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 0.0, 1.0, 1.0))


LINSTC = DEFI_LIST_REEL(VALE=(0.0, 0.5, 1.0))

snl = STAT_NON_LINE(
    CHAM_MATER=AFFMAT,
    MODELE=model,
    EXCIT=(_F(CHARGE=char_cin, FONC_MULT=RAMPE), _F(CHARGE=char_meca, FONC_MULT=RAMPE)),
    INCREMENT=_F(LIST_INST=LINSTC),
    SOLVEUR=_F(METHODE="PETSC", PRE_COND="SANS", RESI_RELA=1.0e-10),
)

TEST_RESU(
    RESU=(
        _F(
            GROUP_NO=("N1",),
            NOM_CHAM="DEPL",
            NOM_CMP="DX",
            INST=1.0,
            RESULTAT=snl,
            VALE_CALC=0.979591885789513,
        ),
        _F(
            GROUP_NO=("N3",),
            NOM_CHAM="DEPL",
            NOM_CMP="DX",
            INST=1.0,
            RESULTAT=snl,
            VALE_CALC=1.020408114214617,
        ),
    )
)

FIN()
