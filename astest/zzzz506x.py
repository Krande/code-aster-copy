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
from code_aster.Utilities import haveMPI

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"), INFO=1)

test = CA.TestCase()

mesh = LIRE_MAILLAGE(FORMAT="MED", UNITE=20, PARTITIONNEUR="PTSCOTCH")

model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

FTRACTUB = DEFI_FONCTION(
    NOM_PARA="EPSI",
    VALE=(1000.0, 2000.0, 2000.0, 5000.0),
    PROL_DROITE="LINEAIRE",
    PROL_GAUCHE="LINEAIRE",
)

# Very high elasticity limit to simulate elasticity
plasticity = CREA_LIB_MFRONT(NOM_COMPOR="Plasticity", UNITE_MFRONT=38)

acier = DEFI_MATERIAU(
    ELAS=_F(E=200000.0, NU=0.3),
    ECRO_LINE=_F(D_SIGM_EPSI=2000.0, SY=200000.0),
    MFRONT=_F(LISTE_COEF=(200000.0, 2000.0, 200000.0, 0.3)),
)

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=acier))


encast = AFFE_CHAR_CINE(MODELE=model, MECA_IMPO=(_F(GROUP_MA="BAS", DX=0, DY=0.0, DZ=0.0),))

depl = AFFE_CHAR_CINE(MODELE=model, MECA_IMPO=(_F(GROUP_MA="HAUT", DZ=1.0),))

LIST = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=2))

RAMPE = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 0.0, 1000.0, 1000.0))

if haveMPI():
    linear_solver = {"METHODE": "PETSC", "PRE_COND": "GAMG", "RESI_RELA": 1e-13}
else:
    linear_solver = {"METHODE": "MUMPS", "RENUM": "AUTO", "NPREC": 8}

# STAT_NON_LINE DE REFERENCE
SOLUT = STAT_NON_LINE(
    MODELE=model,
    CHAM_MATER=mater,
    EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE), _F(CHARGE=depl, FONC_MULT=RAMPE)),
    COMPORTEMENT=_F(RELATION="MFRONT", COMPOR_MFRONT=plasticity),
    CONVERGENCE=_F(RESI_GLOB_MAXI=5e-9),
    NEWTON=_F(REAC_INCR=1, PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    INCREMENT=_F(LIST_INST=LIST),
    SOLVEUR=linear_solver,
    INFO=1,
)

print("Field in original SNL:", flush=True)
SOLUT.printListOfFields()

# NEW STAT_NON_LINE
SOLUN = MECA_NON_LINE(
    MODELE=model,
    CHAM_MATER=mater,
    EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE), _F(CHARGE=depl, FONC_MULT=RAMPE)),
    COMPORTEMENT=_F(RELATION="MFRONT", COMPOR_MFRONT=plasticity),
    CONVERGENCE=_F(RESI_GLOB_MAXI=5e-9),
    NEWTON=_F(REAC_INCR=1, PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    INCREMENT=_F(LIST_INST=LIST),
    SOLVEUR=linear_solver,
    INFO=1,
)

print("Field in new SNL:", flush=True)
SOLUN.printListOfFields()

# =======================================================
#             IMPR_RESU et LIRE_RESU
# =======================================================

# =========================================================
#          DETERMINATION DE LA REFERENCE
# =========================================================

nbIndexes = SOLUT.getNumberOfIndexes()
test.assertEqual(SOLUT.getNumberOfIndexes(), SOLUN.getNumberOfIndexes())


# =========================================================
#            REALISATION DES TESTS
# =========================================================

for rank in range(nbIndexes):
    # ON EXTRAIT LES CHAMPS A TESTER au dernier instant
    DEPL_REF = SOLUT.getField("DEPL", rank)
    SIGMA_REF = SOLUT.getField("SIEF_ELGA", rank)
    VARI_REF = SOLUT.getField("VARI_ELGA", rank)

    DEPL = SOLUN.getField("DEPL", rank)
    SIGMA = SOLUN.getField("SIEF_ELGA", rank)
    VARI = SOLUN.getField("VARI_ELGA", rank)

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


FIN()
