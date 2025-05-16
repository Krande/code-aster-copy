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
from code_aster.Utilities import config

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"), INFO=1)

mesh = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

acier = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3), ECRO_LINE=_F(D_SIGM_EPSI=2000.0, SY=200.0))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=acier))


encast = AFFE_CHAR_MECA(MODELE=model, DDL_IMPO=(_F(GROUP_MA="BAS", DX=0, DY=0.0, DZ=0.0),))

depl = AFFE_CHAR_MECA(MODELE=model, DDL_IMPO=(_F(GROUP_MA="HAUT", DZ=1.0),))

LIST = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=2))

RAMPE = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 0.0, 1000.0, 1000.0))

common_keywords = _F(
    MODELE=model,
    CHAM_MATER=mater,
    EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE), _F(CHARGE=depl, FONC_MULT=RAMPE)),
    COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE", DEFORMATION="GDEF_LOG"),
    NEWTON=_F(REAC_INCR=1, PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    CONVERGENCE=_F(RESI_GLOB_RELA=1e-8),
    INCREMENT=_F(LIST_INST=LIST),
    INFO=1,
)

init_keywords = common_keywords.copy()
init_keywords["INCREMENT"] = _F(LIST_INST=LIST, INST_FIN=0.5)

# STAT_NON_LINE DE REFERENCE
ressnl = STAT_NON_LINE(**init_keywords)
ressnl = STAT_NON_LINE(reuse=ressnl, ETAT_INIT=_F(EVOL_NOLI=ressnl), **common_keywords)

# =========================================================
#          DETERMINATION DE LA REFERENCE
# =========================================================

# ON EXTRAIT LES CHAMPS A TESTER au dernier instant
SIGMA_REF = CREA_CHAMP(
    OPERATION="EXTR", TYPE_CHAM="ELGA_SIEF_R", NOM_CHAM="SIEF_ELGA", RESULTAT=ressnl, INST=1.0
)

VARI_REF = CREA_CHAMP(
    OPERATION="EXTR", TYPE_CHAM="ELGA_VARI_R", NOM_CHAM="VARI_ELGA", RESULTAT=ressnl, INST=1.0
)

TEST_RESU(
    CHAM_ELEM=(
        _F(CRITERE="RELATIF", TYPE_TEST="MIN", CHAM_GD=SIGMA_REF, VALE_CALC=-171.82923889207268),
        _F(CRITERE="RELATIF", TYPE_TEST="MAX", CHAM_GD=SIGMA_REF, VALE_CALC=502.4476860211986),
        _F(CRITERE="RELATIF", TYPE_TEST="MIN", CHAM_GD=VARI_REF, VALE_CALC=-171.76026986621378),
        _F(CRITERE="RELATIF", TYPE_TEST="MAX", CHAM_GD=VARI_REF, VALE_CALC=503.966284799366),
    )
)

# =========================================================
#            SOLVEUR NON LINEAIRE SNES
# =========================================================

if config["ASTER_HAVE_PETSC4PY"]:
    SOLU2 = MECA_NON_LINE(
        MODELE=model,
        CHAM_MATER=mater,
        EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE), _F(CHARGE=depl, FONC_MULT=RAMPE)),
        COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE", DEFORMATION="GDEF_LOG"),
        NEWTON=_F(REAC_INCR=1, PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
        METHODE="NEWTON",
        CONVERGENCE=_F(RESI_GLOB_RELA=1e-8),
        INCREMENT=_F(LIST_INST=LIST),
        SOLVEUR=_F(METHODE="PETSC"),
        INFO=1,
    )

    SIGMA2 = CREA_CHAMP(
        OPERATION="EXTR", TYPE_CHAM="ELGA_SIEF_R", NOM_CHAM="SIEF_ELGA", RESULTAT=SOLU2, INST=1.0
    )

    VARI2 = CREA_CHAMP(
        OPERATION="EXTR", TYPE_CHAM="ELGA_VARI_R", NOM_CHAM="VARI_ELGA", RESULTAT=SOLU2, INST=1.0
    )

    DIF_SIG2 = SIGMA_REF - SIGMA2
    DIF_VAR2 = VARI_REF - VARI2

    TEST_RESU(
        CHAM_ELEM=(
            _F(
                CRITERE="ABSOLU",
                REFERENCE="AUTRE_ASTER",
                PRECISION=1.0e-07,
                TYPE_TEST="MIN",
                CHAM_GD=DIF_SIG2,
                VALE_CALC=1.019122919387e-08,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
            _F(
                CRITERE="ABSOLU",
                REFERENCE="AUTRE_ASTER",
                PRECISION=1.0e-07,
                TYPE_TEST="MAX",
                CHAM_GD=DIF_SIG2,
                VALE_CALC=8.803851869743e-09,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
            _F(
                CRITERE="ABSOLU",
                REFERENCE="AUTRE_ASTER",
                PRECISION=1.0e-07,
                TYPE_TEST="MIN",
                CHAM_GD=DIF_VAR2,
                VALE_CALC=1.018500483951e-08,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
            _F(
                CRITERE="ABSOLU",
                REFERENCE="AUTRE_ASTER",
                PRECISION=1.0e-07,
                TYPE_TEST="MAX",
                CHAM_GD=DIF_VAR2,
                VALE_CALC=8.824088126857e-09,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
        )
    )


FIN()
