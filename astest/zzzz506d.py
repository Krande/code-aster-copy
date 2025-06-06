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

import logging

from code_aster.Commands import *
from code_aster import CA
from code_aster.Utilities import logger

logger.setLevel(logging.DEBUG)

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"), INFO=1)

test = CA.TestCase()

mesh = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

FTRACTUB = DEFI_FONCTION(
    NOM_PARA="EPSI",
    VALE=(1000.0, 2000.0, 2000.0, 5000.0),
    PROL_DROITE="LINEAIRE",
    PROL_GAUCHE="LINEAIRE",
)

# Very high elasticity limit to simulate elasticity
acier = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3), ECRO_LINE=_F(D_SIGM_EPSI=2000.0, SY=200000.0))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=acier))


encast = AFFE_CHAR_MECA(MODELE=model, DDL_IMPO=(_F(GROUP_MA="BAS", DX=0, DY=0.0, DZ=0.0),))

depl = AFFE_CHAR_MECA(MODELE=model, DDL_IMPO=(_F(GROUP_MA="HAUT", DZ=1.0),))

LIST = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=2))

RAMPE = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 0.0, 1000.0, 1000.0))

# STAT_NON_LINE DE REFERENCE
SOLUT = STAT_NON_LINE(
    MODELE=model,
    CHAM_MATER=mater,
    EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE), _F(CHARGE=depl, FONC_MULT=RAMPE)),
    COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE"),
    CONVERGENCE=_F(RESI_GLOB_MAXI=1e-8),
    NEWTON=_F(REAC_INCR=1, PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    INCREMENT=_F(LIST_INST=LIST),
    INFO=1,
)

print("Field in original SNL:", flush=True)
SOLUT.printListOfFields()

# NEW STAT_NON_LINE
SOLUN = MECA_NON_LINE(
    MODELE=model,
    CHAM_MATER=mater,
    EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE), _F(CHARGE=depl, FONC_MULT=RAMPE)),
    COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE"),
    CONVERGENCE=_F(RESI_GLOB_MAXI=1e-8),
    NEWTON=_F(REAC_INCR=1, PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    INCREMENT=_F(LIST_INST=LIST),
    ARCHIVAGE=_F(CHAM_EXCLU=[]),
    INFO=1,
)

print("Field in new SNL:", flush=True)
SOLUN.printListOfFields()

# =======================================================
#             IMPR_RESU et LIRE_RESU
# =======================================================

IMPR_RESU(FORMAT="MED", UNITE=81, RESU=_F(RESULTAT=SOLUT))

IMPR_RESU(FORMAT="MED", UNITE=80, RESU=_F(RESULTAT=SOLUN))

SOLUT = LIRE_RESU(
    MODELE=model,
    FORMAT="MED",
    UNITE=81,
    TYPE_RESU="EVOL_NOLI",
    COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE"),
    CHAM_MATER=mater,
    EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE), _F(CHARGE=depl, FONC_MULT=RAMPE)),
    FORMAT_MED=(
        _F(NOM_RESU="SOLUT", NOM_CHAM="DEPL"),
        _F(NOM_RESU="SOLUT", NOM_CHAM="SIEF_ELGA"),
        _F(NOM_RESU="SOLUT", NOM_CHAM="VARI_ELGA"),
    ),
    TOUT_ORDRE="OUI",
)

SOLUR = LIRE_RESU(
    MODELE=model,
    FORMAT="MED",
    UNITE=80,
    TYPE_RESU="EVOL_NOLI",
    COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE"),
    CHAM_MATER=mater,
    EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE), _F(CHARGE=depl, FONC_MULT=RAMPE)),
    FORMAT_MED=(
        _F(NOM_RESU="SOLUN", NOM_CHAM="DEPL"),
        _F(NOM_RESU="SOLUN", NOM_CHAM="SIEF_ELGA"),
        _F(NOM_RESU="SOLUN", NOM_CHAM="VARI_ELGA"),
    ),
    TOUT_ORDRE="OUI",
)

# =========================================================
#          DETERMINATION DE LA REFERENCE
# =========================================================

nbIndexes = SOLUT.getNumberOfIndexes()
test.assertEqual(SOLUT.getNumberOfIndexes(), SOLUN.getNumberOfIndexes())


# =========================================================
#            REALISATION DES TESTS
# =========================================================

for rank in range(nbIndexes):
    DEPL_REF = SOLUT.getField("DEPL", rank)
    SIGMA_REF = SOLUT.getField("SIEF_ELGA", rank)
    VARI_REF = SOLUT.getField("VARI_ELGA", rank)

    DEPL = SOLUR.getField("DEPL", rank)
    SIGMA = SOLUR.getField("SIEF_ELGA", rank)
    VARI = SOLUR.getField("VARI_ELGA", rank)

    DIF_DEPL = DEPL_REF - DEPL

    DIF_SIG = SIGMA_REF - SIGMA

    DIF_VAR = VARI_REF - VARI

    TEST_RESU(
        CHAM_ELEM=(
            _F(
                CRITERE="ABSOLU",
                REFERENCE="AUTRE_ASTER",
                PRECISION=1.0e-08,
                TYPE_TEST="MIN",
                CHAM_GD=DIF_SIG,
                VALE_CALC=1.5916157281026244e-12,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
            _F(
                CRITERE="ABSOLU",
                REFERENCE="AUTRE_ASTER",
                PRECISION=1.0e-08,
                TYPE_TEST="MAX",
                CHAM_GD=DIF_SIG,
                VALE_CALC=1.1368683772161603e-12,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
            _F(
                CRITERE="ABSOLU",
                REFERENCE="AUTRE_ASTER",
                ORDRE_GRANDEUR=1.0e-08,
                TYPE_TEST="MIN",
                CHAM_GD=DIF_VAR,
                VALE_CALC=0.0,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
            _F(
                CRITERE="ABSOLU",
                REFERENCE="AUTRE_ASTER",
                ORDRE_GRANDEUR=1.0e-8,
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
                REFERENCE="AUTRE_ASTER",
                PRECISION=1.0e-08,
                TYPE_TEST="MIN",
                CHAM_GD=DIF_DEPL,
                VALE_CALC=6.984919309616089e-10,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
            _F(
                CRITERE="ABSOLU",
                REFERENCE="AUTRE_ASTER",
                PRECISION=1.0e-08,
                TYPE_TEST="MAX",
                CHAM_GD=DIF_DEPL,
                VALE_CALC=9.313225746154785e-10,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
        )
    )

# test issue32923
v0 = SOLUR.getField("DEPL", 0)
vmodel = CA.FieldOnNodesReal(model)
vmodel.setValues(1.0)

vnew = vmodel.copyUsingDescription(v0.getDescription())

test.assertSequenceEqual(v0.getComponents(), vnew.getComponents())
test.assertEqual(v0.getDescription().getNumberOfDOFs(), vnew.size())
test.assertAlmostEqual(vnew.norm("NORM_2"), vmodel.norm("NORM_2"), delta=1e-12)


# =========================================================
#           TEST CHAMPS DE RESIDUS
# =========================================================

RESI_GLOB = SOLUN.getField("RESI_NOEU", 2)


RESI_RELA = SOLUN.getField("RESI_RELA_NOEU", 2)

TEST_RESU(
    CHAM_NO=(
        _F(
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            PRECISION=1.0e-08,
            TYPE_TEST="MAX",
            CHAM_GD=RESI_GLOB,
            VALE_CALC=3.026798367500305e-09,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            PRECISION=1.0e-08,
            TYPE_TEST="MAX",
            CHAM_GD=RESI_RELA,
            VALE_CALC=1.1440455182424691e-15,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
    )
)

FIN()
