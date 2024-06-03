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

import code_aster
from code_aster.Commands import *
from code_aster import CA
import numpy as np
from code_aster.CA import MPI
from code_aster.Utilities import PETSc, petscInitialize


CA.init("--test", "--abort")

petscInitialize("-on_error_abort")

test = CA.TestCase()

MA = LIRE_MAILLAGE(PARTITIONNEUR="PTSCOTCH")

MO = AFFE_MODELE(MAILLAGE=MA, AFFE=(_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"),))

MAT = DEFI_MATERIAU(ELAS=_F(E=1000.0, NU=0.3, RHO=1000, AMOR_ALPHA=0.1, AMOR_BETA=0.1))

CHMAT = AFFE_MATERIAU(MAILLAGE=MA, AFFE=_F(TOUT="OUI", MATER=MAT))

sinus = FORMULE(VALE="sin(4*INST*pi*2.)*exp(-10*INST)", NOM_PARA="INST")

x0 = 0.0
y0 = 0.0
z0 = 0.0

x1 = 0.0
y1 = 0.0
z1 = 0.5

ONDE = AFFE_CHAR_MECA_F(
    MODELE=MO,
    ONDE_PLANE=_F(
        DIRECTION=(0.0, 0.0, 1.0),
        TYPE_ONDE="P",
        FONC_SIGNAL=sinus,
        COOR_SOURCE=(x0, y0, z0),
        COOR_REFLECHI=(x1, y1, z1),
        GROUP_MA=("COTE_H",),
    ),
)
BLOQ = AFFE_CHAR_CINE(MODELE=MO, MECA_IMPO=_F(GROUP_MA=("COTE_B",), DX=0, DY=0, DZ=0))
ONDE = AFFE_CHAR_MECA_F(MODELE=MO, PRES_REP=_F(GROUP_MA=("COTE_H",), PRES=sinus), VERI_NORM="NON")

KEL = CALC_MATR_ELEM(OPTION="RIGI_MECA", MODELE=MO, CHAM_MATER=CHMAT)
MEL = CALC_MATR_ELEM(OPTION="MASS_MECA", MODELE=MO, CHAM_MATER=CHMAT)
CEL = CALC_MATR_ELEM(OPTION="AMOR_MECA", MODELE=MO, CHAM_MATER=CHMAT, RIGI_MECA=KEL, MASS_MECA=MEL)
FELEM = CALC_VECT_ELEM(OPTION="CHAR_MECA", CHARGE=ONDE, CHAM_MATER=CHMAT, INST=0.12)
NUMEDDL = NUME_DDL(MATR_RIGI=KEL)


STIFFNESS = ASSE_MATRICE(MATR_ELEM=KEL, NUME_DDL=NUMEDDL, CHAR_CINE=BLOQ)
DAMPING = ASSE_MATRICE(MATR_ELEM=CEL, NUME_DDL=NUMEDDL, CHAR_CINE=BLOQ)
MASS = ASSE_MATRICE(MATR_ELEM=MEL, NUME_DDL=NUMEDDL, CHAR_CINE=BLOQ)

TEMPCAL = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=10, NOMBRE=100))


REFE = DYNA_VIBRA(
    TYPE_CALCUL="TRAN",
    BASE_CALCUL="PHYS",
    MODELE=MO,
    MATR_MASS=MASS,
    MATR_RIGI=STIFFNESS,
    MATR_AMOR=DAMPING,
    EXCIT=(_F(CHARGE=BLOQ), _F(CHARGE=ONDE, COEF_MULT=1.0)),
    INCREMENT=_F(LIST_INST=TEMPCAL),
    INFO=2,
)

ninst = REFE.getNumberOfIndexes()
fref = REFE.getField("DEPL", ninst - 1)

myOpt = "-pc_type asm -sub_pc_type lu -pc_factor_mat_solver_type mumps -ksp_view"
NEW = DYNA_VIBRA(
    TYPE_CALCUL="TRAN",
    BASE_CALCUL="PHYS",
    MODELE=MO,
    MATR_MASS=MASS,
    MATR_RIGI=STIFFNESS,
    MATR_AMOR=DAMPING,
    EXCIT=(_F(CHARGE=BLOQ), _F(CHARGE=ONDE, COEF_MULT=1.0)),
    INCREMENT=_F(LIST_INST=TEMPCAL),
    SOLVEUR=_F(METHODE="PETSC", PRE_COND="SANS", OPTION_PETSC=myOpt),
    INFO=2,
)

ninst = NEW.getNumberOfIndexes()
fnew = NEW.getField("DEPL", ninst - 1)

# -------------------------------

diff = (fref - fnew).norm("NORM_2") / fref.norm("NORM_2")

print(f"""|fref|={fref.norm("NORM_2")} |fnew|={fnew.norm("NORM_2")} """)
test.assertAlmostEqual(fref.norm("NORM_2"), fnew.norm("NORM_2"), places=5)

FIN()
