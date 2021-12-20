# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
from code_aster.Commands import *
from code_aster import MPI

test = code_aster.TestCase()

code_aster.init("--test")

rank = MPI.COMM_WORLD.Get_rank()

pMesh2 = code_aster.ParallelMesh()
pMesh2.readMedFile("zzzz504b/%d.med"%rank, True)

model = AFFE_MODELE(MAILLAGE = pMesh2,
                    AFFE = _F(MODELISATION = "D_PLAN_INCO_UPG",
                              PHENOMENE = "MECANIQUE",
                              TOUT = "OUI",),)

char_cin = AFFE_CHAR_CINE(MODELE=model,
                          MECA_IMPO=(_F(GROUP_NO="N2",
                                        DX=0.,DY=0.,DZ=0.,),
                                     _F(GROUP_NO="N4",
                                        DX=0.,DY=0.,DZ=0.,),),)

char_meca = AFFE_CHAR_MECA(MODELE=model,
                           LIAISON_DDL=_F(GROUP_NO=("N1", "N3"),
                                          DDL=('PRES','PRES'),
                                          COEF_MULT=(1.0,1.0),
                                          COEF_IMPO=0,),
                           DDL_IMPO=(_F(GROUP_NO="N1",PRES=200000.0),
                                          ))
char_meca.debugPrint(10+rank)

MATER1 = DEFI_MATERIAU(ELAS=_F(E=200000.0,
                               NU=0.499999,),)

AFFMAT = AFFE_MATERIAU(MAILLAGE=pMesh2,
                       AFFE=_F(TOUT='OUI',
                               MATER=MATER1,),)

LINSTC=DEFI_LIST_REEL(VALE=( 1.0, 2.0),)


resu = MECA_STATIQUE(CHAM_MATER=AFFMAT,
                     MODELE=model,
                     LIST_INST = LINSTC, INST_FIN=1.0,
                     EXCIT=(_F(CHARGE=char_cin,),
                            _F(CHARGE=char_meca,),),
                     )

ranks = resu.getRanks()
test.assertEqual(len(ranks), 1)
test.assertAlmostEqual(ranks[0], 1.0)
resu.debugPrint(10+rank)

test.assertEqual("resu", resu.userName)
test.assertFalse(resu.hasElementaryCharacteristics())
test.assertFalse(resu.hasElementaryCharacteristics(1))

resu.printMedFile("test"+str(rank)+".med")
#from shutil import copyfile
#copyfile("test"+str(rank)+".med", "/home/siavelis/test"+str(rank)+".med")

MyFieldOnNodes = resu.getFieldOnNodesReal("DEPL", 1)
sfon = MyFieldOnNodes.exportToSimpleFieldOnNodes()
sfon.debugPrint(10+rank)
sfon.updateValuePointers()


# DX displacement on nodes "N1" and "N3", comparison with sequential results
if rank == 0:
    test.assertAlmostEqual(sfon.getValue(1, 0), 1.14977255749554, 6)
elif rank == 1:
    test.assertAlmostEqual(sfon.getValue(1, 0), 1.14977255749554, 6)

FIN()
