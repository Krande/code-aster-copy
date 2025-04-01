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


from code_aster.Commands import *
from code_aster import CA
from code_aster.CA import MPI
from code_aster.Utilities import SharedTmpdir
import os.path as osp

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()
size = MPI.ASTER_COMM_WORLD.Get_size()

mesh = CA.ParallelMesh()
mesh.readMedFile("fort.20")

ACIER = DEFI_MATERIAU(ELAS=_F(RHO=2400.0, NU=0.0, E=6.7e10))

MATER = AFFE_MATERIAU(
    MAILLAGE=mesh, AFFE=(_F(GROUP_MA="VOLUME", MATER=ACIER), _F(GROUP_MA="POUTRE", MATER=ACIER))
)

STRUC_1 = AFFE_MODELE(
    MAILLAGE=mesh,
    AFFE=(
        _F(GROUP_MA="VOLUME", MODELISATION="3D", PHENOMENE="MECANIQUE"),
        _F(GROUP_MA="FACES", MODELISATION="3D", PHENOMENE="MECANIQUE"),
        _F(GROUP_MA="POUTRE", MODELISATION="POU_D_E", PHENOMENE="MECANIQUE"),
    ),
)

CARA_1 = AFFE_CARA_ELEM(
    MODELE=STRUC_1,
    POUTRE=_F(GROUP_MA="POUTRE", SECTION="RECTANGLE", CARA=("HY", "HZ"), VALE=(0.014, 0.014)),
)

FIXA_1 = AFFE_CHAR_MECA(
    MODELE=STRUC_1,
    DDL_IMPO=(
        _F(GROUP_NO="POINT_A", DX=1.0, DY=0.0, DZ=0.0),
        _F(GROUP_NO="DEPL_IMPO", DX=0.0, DY=0.0, DZ=0.0),
    ),
    LIAISON_ELEM=_F(OPTION="3D_POU", GROUP_MA_1="FACES", GROUP_NO_2="NOE1"),
)

stepping = DEFI_LIST_INST(DEFI_LIST=_F(VALE=(0, 1)))

resu = STAT_NON_LINE(
    MODELE=STRUC_1,
    CHAM_MATER=MATER,
    CARA_ELEM=CARA_1,
    EXCIT=_F(CHARGE=FIXA_1),
    INCREMENT=_F(LIST_INST=stepping),
)

resu = CALC_CHAMP(reuse=resu, RESULTAT=resu, CONTRAINTE=("SIEF_NOEU",))

with SharedTmpdir("zzzz503t_") as tmpdir:
    medfile = osp.join(tmpdir.path, "resu_new.med")
    DEFI_FICHIER(UNITE=87, FICHIER=medfile, TYPE="BINARY")
    IMPR_RESU(
        FICHIER_UNIQUE="OUI", PROC0="NON", FORMAT="MED", UNITE=87, RESU=_F(RESULTAT=resu, INST=1.0)
    )
    DEFI_FICHIER(ACTION="LIBERER", UNITE=87)
    test.assertTrue(True)

myFieldOnNodes = resu.getField("DEPL", 1)
sfon = myFieldOnNodes.toSimpleFieldOnNodes()
if rank == 0:
    test.assertAlmostEqual(sfon[2, 0], 0.05555538259520208)
elif rank == 1:
    test.assertAlmostEqual(sfon[6, 0], 0.6666666627174732)

FIN()
