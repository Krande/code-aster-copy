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

# This testcase is a copy of ssnv160f, but using MECA_NON_LINE.
import numpy as np
from scipy.optimize import fsolve

from code_aster.Commands import *
from code_aster import CA

command = MECA_NON_LINE

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>> Isotropic compression test on a 3D HEXA8 element with the CSSM model
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

DEBUT(CODE=_F(NIV_PUB_WEB="INTERNET"), DEBUG=_F(SDVERI="NON"))

### >>>>>>>>>>>>>
### Read the mesh
### <<<<<<<<<<<<<

MAIL = LIRE_MAILLAGE(FORMAT="ASTER")

### >>>>>>>>>>>>>>>>>
### Model affectation
### <<<<<<<<<<<<<<<<<

MODELE = AFFE_MODELE(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

### >>>>>>>>>>>>>>>>>>>
### Material properties
### <<<<<<<<<<<<<<<<<<<

k = 516.0e6
mu = 238.0e6
rho = 0.1
M = 1.38
pc0 = 100.0e3
bt = 30.0
eta = 0.99
om = 32.0
gammahyp = 2.0e-4
nhyp = 0.78
C = 448.0e3

MATER = DEFI_MATERIAU(
    ELAS=_F(
        E=9.0 * k * mu / (3.0 * k + mu), NU=(3.0 * k - 2.0 * mu) / (2.0 * (3.0 * k + mu)), ALPHA=0.0
    ),
    CSSM=_F(
        BulkModulus=k,
        ShearModulus=mu,
        InitCritPress=pc0,
        CritStateSlope=M,
        IncoPlastIndex=bt,
        HypExponent=nhyp,
        HypDistortion=gammahyp,
        MinCritPress=C,
        ShearModulusRatio=rho,
        IsoHardRatio=eta,
        IsoHardIndex=om,
    ),
)

### >>>>>>>>>>>>>>>>>>>>
### Material affectation
### <<<<<<<<<<<<<<<<<<<<

MATE = AFFE_MATERIAU(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", MATER=MATER))

### >>>>>>>>>>
### Time steps
### <<<<<<<<<<

LI1 = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=(_F(JUSQU_A=1.0, NOMBRE=5), _F(JUSQU_A=2.0, NOMBRE=5)))

maxSteps = 22
DEFLIST = DEFI_LIST_INST(
    METHODE="AUTO",
    DEFI_LIST=_F(LIST_INST=LI1, NB_PAS_MAXI=maxSteps),
    ECHEC=_F(EVENEMENT="ERREUR", ACTION="DECOUPE", SUBD_METHODE="AUTO"),
)

### >>>>>>>>>>>>>>>>>>>
### Boundary conditions
### <<<<<<<<<<<<<<<<<<<

# /!\ p is 50 times more than in ssnv160f + reduced ITER_INTE_MAXI to force the timesteps splitting
p = 500.0e4

PRES = DEFI_FONCTION(NOM_PARA="INST", NOM_RESU="PRESSION", VALE=(0.0, 0.0, 1.0, p, 2.0, 100.0))

CHA_PRES = AFFE_CHAR_MECA_F(
    MODELE=MODELE, PRES_REP=_F(GROUP_MA=("HAUT", "DROITE", "ARRIERE"), PRES=PRES), VERI_NORM="OUI"
)

CHA_DEPL = AFFE_CHAR_MECA(
    MODELE=MODELE,
    DDL_IMPO=(
        _F(GROUP_NO=("NO1", "NO2", "NO4", "NO3"), DZ=0.0),
        _F(GROUP_NO=("NO3", "NO4", "NO7", "NO8"), DY=0.0),
        _F(GROUP_NO=("NO2", "NO4", "NO6", "NO8"), DX=0.0),
    ),
)

### >>>>>>>>
### Solution
### <<<<<<<<

RESU = command(
    MODELE=MODELE,
    CHAM_MATER=MATE,
    EXCIT=(_F(CHARGE=CHA_PRES), _F(CHARGE=CHA_DEPL)),
    COMPORTEMENT=_F(RELATION="CSSM", RESI_INTE=1.0e-14, ITER_INTE_MAXI=20),
    INCREMENT=_F(LIST_INST=DEFLIST),
    NEWTON=_F(MATRICE="TANGENTE"),
    CONVERGENCE=_F(ITER_GLOB_MAXI=10, RESI_GLOB_RELA=1.0e-10),
    SOLVEUR=_F(METHODE="MUMPS"),
)

### >>>>>>>>>>>>>>>
### Post-processing
### <<<<<<<<<<<<<<<

RESU = CALC_CHAMP(
    reuse=RESU,
    RESULTAT=RESU,
    CONTRAINTE="SIEF_NOEU",
    DEFORMATION="EPSI_NOEU",
    VARI_INTERNE="VARI_NOEU",
)

### >>>>>
### Tests
### <<<<<


def func(x):
    return p - 2.0 * pc0 * (np.exp(-bt * x) - eta * np.exp(2.0 * om * x))


x0 = -np.log(p / (2.0 * pc0)) / bt
EPVP = fsolve(func, x0)[0]
EPVE = -p / k

test = CA.TestCase()
test.assertEqual(RESU.getNumberOfIndexes() - 1, maxSteps)

TEST_RESU(
    RESU=_F(
        INST=1.0,
        RESULTAT=RESU,
        REFERENCE="ANALYTIQUE",
        NOM_CHAM="EPSI_NOEU",
        GROUP_NO="NO6",
        NOM_CMP="EPXX",
        VALE_REFE=(EPVE + EPVP) / 3.0,
        VALE_CALC=-0.038995719385912855,
    )
)

TEST_RESU(
    RESU=_F(
        INST=1.0,
        RESULTAT=RESU,
        REFERENCE="ANALYTIQUE",
        NOM_CHAM="VARI_NOEU",
        NOM_CMP="V15",
        GROUP_NO="NO6",
        VALE_REFE=EPVP,
        VALE_CALC=-0.10729723567711803,
    )
)

TEST_RESU(
    RESU=_F(
        NUME_ORDRE=maxSteps,
        RESULTAT=RESU,
        NOM_CHAM="EPSI_NOEU",
        GROUP_NO="NO6",
        NOM_CMP="EPXX",
        VALE_CALC=-0.037703755561623496,
    )
)

TEST_RESU(
    RESU=_F(
        NUME_ORDRE=maxSteps,
        RESULTAT=RESU,
        NOM_CHAM="VARI_NOEU",
        NOM_CMP="V15",
        GROUP_NO="NO6",
        VALE_CALC=-0.10729723567711803,
    )
)


FIN()
