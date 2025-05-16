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

# ***********************************************************************
#
# TITRE: ESSAI DE TRACTION AVEC LA LOI DE RANKINE
# MODIFICATION DU TEST SSNV515C
#
# ***********************************************************************

# ======================================================================
DEBUT(CODE="OUI")

# ***********************************************************************
#
#    MAILLAGE ET MODELE
#
# ***********************************************************************

MAILLAGE = LIRE_MAILLAGE(FORMAT="MED", PARTITIONNEUR="PTSCOTCH")

MODELE = AFFE_MODELE(
    MAILLAGE=MAILLAGE, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="AXIS")
)

MAILLAGE = MODI_MAILLAGE(
    reuse=MAILLAGE,
    MAILLAGE=MAILLAGE,
    ORIE_PEAU=_F(GROUP_MA_PEAU=("DROIT", "GAUCHE", "BAS", "HAUT")),
    INFO=1,
)

# ***********************************************************************
#
#    LISTE D'INSTANTS
#
# ***********************************************************************

tarret = 10.0
tfin = 20.0
npas = 30
temps_max = 30.0

dtemps = temps_max / npas

ltemps = [dtemps * i for i in range(npas + 1)]

TEMPS = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=(_F(JUSQU_A=20, NOMBRE=10),))
TEMPS2 = DEFI_LIST_REEL(DEBUT=20, INTERVALLE=(_F(JUSQU_A=temps_max, NOMBRE=5),))

# ***********************************************************************
#
#    DONNEES MATERIAU
#
# ***********************************************************************

# modules mecaniques [kPa]
# ------------------------
# K=516.2E6
# G=238.2E6
# # =>
# YOUNG = 9.*K*G /(3.*K+G)
# POISSON = (3.*K-2.*G) /(6.*K+2.*G)

YOUNG = 1e6
POISSON = 0.25
SIGMA_T = 1.0e3
K = YOUNG / 3.0 / (1.0 - 2.0 * POISSON)
G = YOUNG / 2.0 / (1.0 + POISSON)

SOL = DEFI_MATERIAU(ELAS=_F(E=YOUNG, NU=POISSON, ALPHA=0.0), RANKINE=_F(SIGMA_T=SIGMA_T), INFO=1)

CHMAT = AFFE_MATERIAU(MAILLAGE=MAILLAGE, AFFE=_F(TOUT="OUI", MATER=SOL))

# ***********************************************************************
#
#    CHARGEMENTS
#
# ***********************************************************************

# pression de preconsolidation [en kPa]
P0 = 5.0e3
EPZZ = 0.03

npas = 300
dtemps = temps_max / npas
linst = [dtemps * i for i in range(npas)]

SIGLAT = AFFE_CHAR_MECA(MODELE=MODELE, PRES_REP=_F(GROUP_MA=("DROIT",), PRES=P0))

DEPHAUT = AFFE_CHAR_CINE(MODELE=MODELE, MECA_IMPO=(_F(GROUP_MA=("HAUT",), DY=1.0),))

DEPL_1 = AFFE_CHAR_CINE(
    MODELE=MODELE, MECA_IMPO=(_F(GROUP_MA="BAS", DY=0.0), _F(GROUP_NO="B", DX=0.0))
)


COEF2 = DEFI_FONCTION(
    NOM_PARA="INST", PROL_DROITE="CONSTANT", VALE=(0.0, 0.0, 20, EPZZ, temps_max, 0)
)

COEF3 = DEFI_FONCTION(NOM_PARA="INST", PROL_DROITE="CONSTANT", VALE=(0.0, 1.0))


# ***********************************************************************
#
#    PRECONSOLIDATION INITIALE A 5KPA
#
# ***********************************************************************

SIG0 = CREA_CHAMP(
    INFO=2,
    TYPE_CHAM="ELGA_SIEF_R",
    OPERATION="AFFE",
    MODELE=MODELE,
    PROL_ZERO="OUI",
    AFFE=_F(GROUP_MA="BLOC", NOM_CMP=("SIXX", "SIYY", "SIZZ"), VALE=(-P0, -P0, -P0)),
)


# ***********************************************************************
#
#    ESSAI TRIAXIAL DRAINE
#
# ***********************************************************************

U1 = STAT_NON_LINE(
    MODELE=MODELE,
    CHAM_MATER=CHMAT,
    EXCIT=(
        _F(CHARGE=SIGLAT, FONC_MULT=COEF3),
        _F(CHARGE=DEPHAUT, FONC_MULT=COEF2),
        _F(CHARGE=DEPL_1, FONC_MULT=COEF3),
    ),
    ETAT_INIT=_F(SIGM=SIG0),
    COMPORTEMENT=_F(RELATION="RANKINE"),
    NEWTON=_F(MATRICE="TANGENTE", PREDICTION="ELASTIQUE", REAC_ITER=1),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=30, ARRET="OUI"),
    SOLVEUR=_F(METHODE="MUMPS", NPREC=8),
    INCREMENT=_F(LIST_INST=TEMPS),
)

U2 = MECA_NON_LINE(
    MODELE=MODELE,
    CHAM_MATER=CHMAT,
    EXCIT=(
        _F(CHARGE=SIGLAT, FONC_MULT=COEF3),
        _F(CHARGE=DEPHAUT, FONC_MULT=COEF2),
        _F(CHARGE=DEPL_1, FONC_MULT=COEF3),
    ),
    ETAT_INIT=_F(SIGM=SIG0),
    COMPORTEMENT=_F(RELATION="RANKINE"),
    NEWTON=_F(MATRICE="TANGENTE", PREDICTION="ELASTIQUE", REAC_ITER=1),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=30, ARRET="OUI"),
    SOLVEUR=_F(METHODE="MUMPS", NPREC=8),
    INCREMENT=_F(LIST_INST=TEMPS),
)

# TEST ETAT_INIT  AVEC SIGM
# *************************

# ON EXTRAIT LES CHAMPS A TESTER au dernier instant

nbIndexes = U1.getNumberOfIndexes()

DEPL_REF = U1.getField("DEPL", nbIndexes - 1)
SIGMA_REF = U1.getField("SIEF_ELGA", nbIndexes - 1)
VARI_REF = U1.getField("VARI_ELGA", nbIndexes - 1)

DEPL = U2.getField("DEPL", nbIndexes - 1)
SIGMA = U2.getField("SIEF_ELGA", nbIndexes - 1)
VARI = U2.getField("VARI_ELGA", nbIndexes - 1)

DIF_DEPL = DEPL_REF - DEPL
DIF_SIG = SIGMA_REF - SIGMA
DIF_VAR = VARI_REF - VARI

TEST_RESU(
    CHAM_ELEM=(
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            TYPE_TEST="MIN",
            CHAM_GD=DIF_SIG,
            VALE_CALC=1e-10,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_SIG,
            VALE_CALC=1e-10,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            ORDRE_GRANDEUR=10e-3,
            TYPE_TEST="MIN",
            CHAM_GD=DIF_VAR,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            ORDRE_GRANDEUR=10e-3,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_VAR,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
    )
)

# EVOL_NOLI
# *********

U3 = STAT_NON_LINE(
    MODELE=MODELE,
    CHAM_MATER=CHMAT,
    EXCIT=(
        _F(CHARGE=SIGLAT, FONC_MULT=COEF3),
        _F(CHARGE=DEPHAUT, FONC_MULT=COEF2),
        _F(CHARGE=DEPL_1, FONC_MULT=COEF3),
    ),
    ETAT_INIT=_F(EVOL_NOLI=U1),
    COMPORTEMENT=_F(RELATION="RANKINE"),
    NEWTON=_F(MATRICE="TANGENTE", PREDICTION="ELASTIQUE", REAC_ITER=1),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=10, ARRET="OUI"),
    SOLVEUR=_F(METHODE="MUMPS", NPREC=8),
    INCREMENT=_F(LIST_INST=TEMPS2, INST_INIT=20.0),
)

U4 = MECA_NON_LINE(
    MODELE=MODELE,
    CHAM_MATER=CHMAT,
    EXCIT=(
        _F(CHARGE=SIGLAT, FONC_MULT=COEF3),
        _F(CHARGE=DEPHAUT, FONC_MULT=COEF2),
        _F(CHARGE=DEPL_1, FONC_MULT=COEF3),
    ),
    ETAT_INIT=_F(EVOL_NOLI=U2),
    COMPORTEMENT=_F(RELATION="RANKINE"),
    NEWTON=_F(MATRICE="TANGENTE", PREDICTION="ELASTIQUE", REAC_ITER=1),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-8, ITER_GLOB_MAXI=10, ARRET="OUI"),
    SOLVEUR=_F(METHODE="MUMPS", NPREC=8),
    INCREMENT=_F(LIST_INST=TEMPS2, INST_INIT=20.0),
)

# TEST EVOL_NOLI
# *******************************

# Test des intervalles de temps
test = CA.TestCase()

nbIndexes3 = U3.getNumberOfIndexes()
nbIndexes4 = U4.getNumberOfIndexes()

range3 = [U3.getTime(i) for i in range(nbIndexes3)]
range4 = [U3.getTime(i) for i in range(nbIndexes4)]

test.assertEqual(nbIndexes3, nbIndexes4)
test.assertEqual(range3, range4)

DEPL_REF = U3.getField("DEPL", nbIndexes3 - 1)
SIGMA_REF = U3.getField("SIEF_ELGA", nbIndexes3 - 1)
VARI_REF = U3.getField("VARI_ELGA", nbIndexes3 - 1)

DEPL = U4.getField("DEPL", nbIndexes3 - 1)
SIGMA = U4.getField("SIEF_ELGA", nbIndexes3 - 1)
VARI = U4.getField("VARI_ELGA", nbIndexes3 - 1)

DIF_DEPL = DEPL_REF - DEPL
DIF_SIG = SIGMA_REF - SIGMA
DIF_VAR = VARI_REF - VARI

TEST_RESU(
    CHAM_ELEM=(
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            TYPE_TEST="MIN",
            CHAM_GD=DIF_SIG,
            VALE_CALC=1e-10,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_SIG,
            VALE_CALC=1e-10,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            ORDRE_GRANDEUR=10e-3,
            TYPE_TEST="MIN",
            CHAM_GD=DIF_VAR,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            ORDRE_GRANDEUR=10e-3,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_VAR,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
    )
)

TEST_RESU(
    CHAM_NO=(
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-8,
            TYPE_TEST="MIN",
            CHAM_GD=DIF_DEPL,
            VALE_CALC=1.0e-10,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-8,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_DEPL,
            VALE_CALC=1.0e-10,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
    )
)

#### TEST ETAT_INIT U 0

U0 = CREA_CHAMP(
    INFO=2,
    TYPE_CHAM="NOEU_DEPL_R",
    OPERATION="AFFE",
    MODELE=MODELE,
    AFFE=_F(GROUP_MA="BLOC", NOM_CMP=("DX", "DY", "DZ"), VALE=(-0.001, -0.001, -0.001)),
)

U5 = STAT_NON_LINE(
    MODELE=MODELE,
    CHAM_MATER=CHMAT,
    EXCIT=(
        _F(CHARGE=SIGLAT, FONC_MULT=COEF3),
        _F(CHARGE=DEPHAUT, FONC_MULT=COEF2),
        _F(CHARGE=DEPL_1, FONC_MULT=COEF3),
    ),
    ETAT_INIT=_F(SIGM=SIG0, DEPL=U0),
    COMPORTEMENT=_F(RELATION="RANKINE"),
    NEWTON=_F(MATRICE="TANGENTE", PREDICTION="ELASTIQUE", REAC_ITER=1),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=30, ARRET="OUI"),
    SOLVEUR=_F(METHODE="MUMPS", NPREC=8),
    INCREMENT=_F(LIST_INST=TEMPS, INST_FIN=2.0),
)

U6 = MECA_NON_LINE(
    MODELE=MODELE,
    CHAM_MATER=CHMAT,
    EXCIT=(
        _F(CHARGE=SIGLAT, FONC_MULT=COEF3),
        _F(CHARGE=DEPHAUT, FONC_MULT=COEF2),
        _F(CHARGE=DEPL_1, FONC_MULT=COEF3),
    ),
    ETAT_INIT=_F(SIGM=SIG0, DEPL=U0),
    COMPORTEMENT=_F(RELATION="RANKINE"),
    NEWTON=_F(MATRICE="TANGENTE", PREDICTION="ELASTIQUE", REAC_ITER=1),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=30, ARRET="OUI"),
    SOLVEUR=_F(METHODE="MUMPS", NPREC=8),
    INCREMENT=_F(LIST_INST=TEMPS, INST_FIN=2.0),
)

SIGMA_REF = U5.getField("SIEF_ELGA", 1)
SIGMA = U6.getField("SIEF_ELGA", 1)

DIF_SIG = SIGMA_REF - SIGMA

TEST_RESU(
    CHAM_ELEM=(
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            TYPE_TEST="MIN",
            CHAM_GD=DIF_SIG,
            VALE_CALC=1e-10,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
    )
)

FIN()
