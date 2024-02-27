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
pMesh2.readMedFile("zzzz504h/%d.med" % rank, partitioned=True)

model = AFFE_MODELE(
    MAILLAGE=pMesh2, AFFE=_F(MODELISATION="POU_D_T", PHENOMENE="MECANIQUE", TOUT="OUI")
)

cara_elem = AFFE_CARA_ELEM(
    MODELE=model, POUTRE=_F(GROUP_MA="poutre", SECTION="CERCLE", CARA="R", VALE=0.3)
)

char_meca = AFFE_CHAR_MECA(
    MODELE=model, DDL_POUTRE=_F(ANGL_VRIL=0.0, DX=0.0, DY=0.0, DZ=0.0, GROUP_NO="N1")
)
char_meca2 = AFFE_CHAR_MECA(MODELE=model, FORCE_NODALE=_F(GROUP_NO="N3", FX=10000.0))
char_meca.debugPrint(10 + rank)

MATER1 = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3))

AFFMAT = AFFE_MATERIAU(MAILLAGE=pMesh2, AFFE=_F(TOUT="OUI", MATER=MATER1))

resu = MECA_STATIQUE(
    CHAM_MATER=AFFMAT,
    MODELE=model,
    CARA_ELEM=cara_elem,
    EXCIT=(_F(CHARGE=char_meca), _F(CHARGE=char_meca2)),
    SOLVEUR=_F(METHODE="PETSC", PRE_COND="SANS", RESI_RELA=1.0e-10),
)
resu.debugPrint(10 + rank)

test.assertTrue(resu.hasElementaryCharacteristics())
test.assertTrue(resu.hasElementaryCharacteristics(1))
test.assertFalse(resu.hasElementaryCharacteristics(4))

resu.printMedFile("test" + str(rank) + ".med")
# from shutil import copyfile
# copyfile("test"+str(rank)+".med", "/home/siavelis/test"+str(rank)+".med")

MyFieldOnNodes = resu.getField("DEPL", 1)
sfon = MyFieldOnNodes.toSimpleFieldOnNodes()
sfon.debugPrint(10 + rank)
sfon.build()

# DX displacement on nodes "N1" and "N3", comparison with sequential results
if rank == 0:
    test.assertAlmostEqual(sfon[0, 0], 0.0)
elif rank == 1:
    test.assertAlmostEqual(sfon[0, 0], 0.5305164769729844)

FIN()
