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
pMesh2.readMedFile("zzzz504k/%d.med" % rank, partitioned=True)

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
    LIAISON_DDL=_F(GROUP_NO=("N3", "N5"), DDL=("DX", "DY"), COEF_MULT=(1.0, -1.0), COEF_IMPO=0),
    DDL_IMPO=_F(GROUP_NO="N3", DX=1.0),
)
char_meca.debugPrint(10 + rank)

MATER1 = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3))

AFFMAT = AFFE_MATERIAU(MAILLAGE=pMesh2, AFFE=_F(TOUT="OUI", MATER=MATER1))

resu = MECA_STATIQUE(
    CHAM_MATER=AFFMAT,
    MODELE=model,
    EXCIT=(_F(CHARGE=char_cin), _F(CHARGE=char_meca)),
    SOLVEUR=_F(METHODE="PETSC", PRE_COND="SANS"),
)
resu.debugPrint(10 + rank)

resu.printMedFile("test" + str(rank) + ".med")

MyFieldOnNodes = resu.getField("DEPL", 1)
sfon = MyFieldOnNodes.toSimpleFieldOnNodes()
sfon.debugPrint(10 + rank)
sfon.build()

resu1 = [0.6980073863837113, 1.0]
test.assertAlmostEqual(sfon[1, 0], resu1[rank])
resu2 = [0.5461603079530792, 0.5461603079530792]
test.assertAlmostEqual(sfon[2, 0], resu2[rank])

FIN()
