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

DEBUT(CODE="OUI", IGNORE_ALARM="MODELISA8_14")

from code_aster.CA import MPI, TestCase

test = TestCase()

MA = LIRE_MAILLAGE(FORMAT="MED", PARTITIONNEUR="PTSCOTCH")

MOD = AFFE_MODELE(
    MAILLAGE=MA,
    AFFE=(
        _F(GROUP_MA=("Beton", "Encast", "Press"), PHENOMENE="MECANIQUE", MODELISATION="3D"),
        _F(GROUP_MA=(f"Cable{0}", f"Cable{1}"), PHENOMENE="MECANIQUE", MODELISATION="BARRE"),
    ),
)

# Orientation de tous les elements de surface
MA = MODI_MAILLAGE(reuse=MA, MAILLAGE=MA, ORIE_PEAU_3D=_F(GROUP_MA="Press"))

# Definition et effectation des materiaux

BTN_GEN = DEFI_MATERIAU(
    ELAS=_F(E=30000e6, NU=0.2, RHO=2500.0), BPEL_BETON=_F(PERT_FLUA=0.0, PERT_RETR=0.0)
)

ACI_CAB = DEFI_MATERIAU(
    ELAS=_F(E=200000e6, NU=0.0, RHO=7800), BPEL_ACIER=_F(RELAX_1000=0.0, F_PRG=0.0)
)


MATER = AFFE_MATERIAU(
    MAILLAGE=MA,
    AFFE=(
        _F(GROUP_MA="Beton", MATER=BTN_GEN),
        _F(GROUP_MA=(f"Cable{0}", f"Cable{1}"), MATER=ACI_CAB),
    ),
)
# Definition et affectation des caracteristiques des elements de structure
ELEM = AFFE_CARA_ELEM(
    MODELE=MOD,
    BARRE=_F(GROUP_MA=(f"Cable{0}", f"Cable{1}"), SECTION="CERCLE", CARA="R", VALE=0.005),
)

LIAISON_0 = AFFE_CHAR_MECA(
    MODELE=MOD,
    DOUBLE_LAGRANGE="NON",
    LIAISON_MAIL=_F(GROUP_MA_MAIT=f"Cable{0}_vol", GROUP_MA_ESCL=f"Cable{0}"),
)

LIAISON_1 = AFFE_CHAR_MECA(
    MODELE=MOD,
    DOUBLE_LAGRANGE="NON",
    LIAISON_MAIL=_F(GROUP_MA_MAIT=f"Cable{1}_vol", GROUP_MA_ESCL=f"Cable{1}"),
)


GRAV = AFFE_CHAR_MECA(
    MODELE=MOD, DOUBLE_LAGRANGE="NON", PESANTEUR=_F(GRAVITE=9.81, DIRECTION=(0.0, 0.0, -1.0))
)

APP = AFFE_CHAR_CINE(MODELE=MOD, MECA_IMPO=_F(GROUP_MA="Encast", DX=0.0, DY=0.0, DZ=0.0))

PRESS = AFFE_CHAR_MECA(
    MODELE=MOD, DOUBLE_LAGRANGE="NON", PRES_REP=(_F(GROUP_MA="Press", PRES=1e6),)
)

RESU = MECA_STATIQUE(
    MODELE=MOD,
    CHAM_MATER=MATER,
    CARA_ELEM=ELEM,
    EXCIT=(_F(CHARGE=APP), _F(CHARGE=PRESS), _F(CHARGE=LIAISON_0), _F(CHARGE=LIAISON_1)),
)

depl = RESU.getField("DEPL", 1).getValues()
d_min = MPI.ASTER_COMM_WORLD.allreduce(min(depl), MPI.MIN)
d_max = MPI.ASTER_COMM_WORLD.allreduce(max(depl), MPI.MAX)
# print("depl max :", d_max)
# print("depl min :", d_min)
test.assertAlmostEqual(d_min, -0.0005074289446608608)
test.assertAlmostEqual(d_max, 2.279187067594351e-06)

FIN()
