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


CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

if MPI.ASTER_COMM_WORLD.Get_size() > 1:
    is_parallel = True
else:
    is_parallel = False

if is_parallel:
    MA = CA.ParallelMesh()
else:
    MA = CA.Mesh()

MA.readMedFile("zzzz502s.mmed")

MO = AFFE_MODELE(
    MAILLAGE=MA,
    DISTRIBUTION=_F(METHODE="CENTRALISE"),
    AFFE=(
        _F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"),
        _F(GROUP_MA=("COTE_H"), PHENOMENE="MECANIQUE", MODELISATION="3D_ABSO"),
    ),
)

MAT = DEFI_MATERIAU(ELAS=_F(E=1000.0, NU=0.3, RHO=1000, AMOR_ALPHA=0.1, AMOR_BETA=0.1))

CHMAT = AFFE_MATERIAU(MAILLAGE=MA, AFFE=_F(TOUT="OUI", MATER=MAT))

sinus = FORMULE(VALE="sin(4*INST*pi*2.)", NOM_PARA="INST")

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

KEL = CALC_MATR_ELEM(OPTION="RIGI_MECA", MODELE=MO, CHAM_MATER=CHMAT)
MEL = CALC_MATR_ELEM(OPTION="MASS_MECA", MODELE=MO, CHAM_MATER=CHMAT)
CEL = CALC_MATR_ELEM(OPTION="AMOR_MECA", MODELE=MO, CHAM_MATER=CHMAT, RIGI_MECA=KEL, MASS_MECA=MEL)


NUMEDDL = NUME_DDL(MATR_RIGI=KEL)


STIFFNESS = ASSE_MATRICE(MATR_ELEM=KEL, NUME_DDL=NUMEDDL)

DAMPING = ASSE_MATRICE(MATR_ELEM=CEL, NUME_DDL=NUMEDDL)

MASS = ASSE_MATRICE(MATR_ELEM=MEL, NUME_DDL=NUMEDDL)

LISTINST = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=(_F(JUSQU_A=1, NOMBRE=10),))

# Algebraic Multigrid Preconditioner
DYNA = DYNA_VIBRA(
    TYPE_CALCUL="TRAN",
    BASE_CALCUL="PHYS",
    MODELE=MO,
    CHAM_MATER=CHMAT,
    MATR_MASS=MASS,
    MATR_RIGI=STIFFNESS,
    MATR_AMOR=DAMPING,
    EXCIT=(_F(CHARGE=ONDE),),
    OBSERVATION=(
        _F(
            CRITERE="RELATIF",
            EVAL_CHAM="VALE",
            EVAL_CMP="VALE",
            GROUP_MA="COTE_H",
            NOM_CHAM="DEPL",
            NOM_CMP="DX",
            OBSE_ETAT_INIT="OUI",
            PAS_OBSE=1,
            PRECISION=1e-06,
        ),
        _F(
            CRITERE="RELATIF",
            EVAL_CHAM="VALE",
            EVAL_CMP="VALE",
            GROUP_MA="COTE_B",
            NOM_CHAM="DEPL",
            NOM_CMP="DX",
            OBSE_ETAT_INIT="OUI",
            PAS_OBSE=1,
            PRECISION=1e-06,
        ),
    ),
    SOLVEUR=_F(METHODE="PETSC", PRE_COND="GAMG", RESI_RELA=1e-10),
    INCREMENT=_F(LIST_INST=LISTINST),
    ARCHIVAGE=_F(LIST_INST=LISTINST),
    SCHEMA_TEMPS=_F(SCHEMA="NEWMARK"),
    INFO=2,
)

DEPL = DYNA.getField("DEPL", 2)
DEPL.setMesh(MA)

# sequential comparaison
seq_value = 15.652508028144187
mpi_value = DEPL.norm("NORM_1")
test.assertTrue(abs(seq_value - mpi_value) / abs(seq_value) < 1e-6)


# Domain Decomposition Preconditioner
DYNA = DYNA_VIBRA(
    TYPE_CALCUL="TRAN",
    BASE_CALCUL="PHYS",
    MODELE=MO,
    CHAM_MATER=CHMAT,
    MATR_MASS=MASS,
    MATR_RIGI=STIFFNESS,
    MATR_AMOR=DAMPING,
    EXCIT=(_F(CHARGE=ONDE),),
    OBSERVATION=(
        _F(
            CRITERE="RELATIF",
            EVAL_CHAM="VALE",
            EVAL_CMP="VALE",
            GROUP_MA="COTE_H",
            NOM_CHAM="DEPL",
            NOM_CMP="DX",
            OBSE_ETAT_INIT="OUI",
            PAS_OBSE=1,
            PRECISION=1e-06,
        ),
        _F(
            CRITERE="RELATIF",
            EVAL_CHAM="VALE",
            EVAL_CMP="VALE",
            GROUP_MA="COTE_B",
            NOM_CHAM="DEPL",
            NOM_CMP="DX",
            OBSE_ETAT_INIT="OUI",
            PAS_OBSE=1,
            PRECISION=1e-06,
        ),
    ),
    SOLVEUR=_F(METHODE="PETSC", PRE_COND="HPDDM", RESI_RELA=1e-10),
    INCREMENT=_F(LIST_INST=LISTINST),
    ARCHIVAGE=_F(LIST_INST=LISTINST),
    SCHEMA_TEMPS=_F(SCHEMA="NEWMARK"),
    INFO=2,
)

DEPL = DYNA.getField("DEPL", 2)
DEPL.setMesh(MA)

# sequential comparaison
seq_value = 15.652508028144187
mpi_value = DEPL.norm("NORM_1")
test.assertTrue(abs(seq_value - mpi_value) / abs(seq_value) < 1e-6)


### On convertit le chargement onde_plane en transitoire dynamique temporel

CH_ON = CREA_RESU(
    OPERATION="CONV_CHAR",
    TYPE_RESU="DYNA_TRANS",
    CONV_CHAR=_F(
        CHAM_MATER=CHMAT,
        MATR_RIGI=STIFFNESS,
        CHARGE=(ONDE,),
        PRECISION=1.0e-6,
        CRITERE="RELATIF",
        LIST_INST=LISTINST,
    ),
)


### Passage temporel-fréquentiel pour CL
CHAONF = REST_SPEC_TEMP(
    RESULTAT=CH_ON,
    METHODE="PROL_ZERO",
    SYMETRIE="NON",
    NOM_CHAM="DEPL",
    N_PUIS=0,
    INFO=2,
    ACCELERATION_MPI="OUI",
)

ref = [
    [1.581553270295848e-12, 0.0],
    [167.95544903379576, 69.56942486426169],
    [112.4530874891084, 112.45308748910587],
    [41.9033036594456, 101.16352400287337],
    [421.4813799041548, 2.0106544693134724e-12],
    [117.12937779837777, 282.77533243316634],
    [1672.9706890353066, 1672.9706890353102],
    [2533.8932869201917, 1049.5729650484825],
]

for i, idx in enumerate(CHAONF.getIndexes()):
    f = CHAONF.getField("DEPL", idx)
    test.assertAlmostEqual(f.getRealPart().norm("NORM_1"), ref[i][0])
    test.assertAlmostEqual(f.getImaginaryPart().norm("NORM_1"), ref[i][1])


FIN()
