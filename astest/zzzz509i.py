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

test = CA.TestCase()

DEBUT(CODE="OUI")

# mesh and properties extracted from forma03
mesh = LIRE_MAILLAGE(FORMAT="MED")

mesh = MODI_MAILLAGE(reuse=mesh, MAILLAGE=mesh, ORIE_PEAU=_F(GROUP_MA_PEAU="haut"))

model = AFFE_MODELE(
    MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="C_PLAN")
)

curve = LIRE_FONCTION(UNITE=21, NOM_PARA="EPSI", PROL_DROITE="CONSTANT")

steel = DEFI_MATERIAU(
    ELAS=_F(E=200000.0, NU=0.3), TRACTION=_F(SIGM=curve), ECRO_LINE=_F(D_SIGM_EPSI=1930.0, SY=181.0)
)

material = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=steel))

sym_bottom = AFFE_CHAR_CINE(MODELE=model, MECA_IMPO=(_F(GROUP_MA="bas", DY=0.0)))
sym_left = AFFE_CHAR_CINE(MODELE=model, MECA_IMPO=_F(GROUP_MA="gauche", DX=0.0))

load = AFFE_CHAR_MECA(MODELE=model, FORCE_CONTOUR=_F(GROUP_MA="haut", FY=1.0))

mult_func = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 0.0, 200.0, 1000.0))

time_start = 0.0
time_inter = 20.0
time_end = 30.0

time_values = DEFI_LIST_REEL(
    DEBUT=time_start, INTERVALLE=(_F(JUSQU_A=time_inter, NOMBRE=4), _F(JUSQU_A=time_end, NOMBRE=8))
)
times = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=time_values), ECHEC=_F(SUBD_NIVEAU=3))
# with >=50, it raises FACTOR_60...

# with test.assertRaisesRegex(CA.AsterError, "niveaux.*subdivision.*atteint"):
with test.assertRaisesRegex(
    CA.AsterError, "Le nombre maximal autorisé SUBD_NIVEAU.*de niveaux de subdivision est dépassé."
):
    MECA_NON_LINE(
        INCREMENT=_F(LIST_INST=times, INST_FIN=time_inter),
        MODELE=model,
        CHAM_MATER=material,
        EXCIT=(_F(CHARGE=sym_bottom), _F(CHARGE=load, FONC_MULT=mult_func)),
        COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE"),
    )

test.printSummary()

FIN()
