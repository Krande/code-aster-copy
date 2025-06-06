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

# This testcase is a copy of ssnp15f, but using MECA_NON_LINE.

from code_aster.Commands import *

command = MECA_NON_LINE

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="NON"))

# ......................................................................
# MODELISATION F : IDEM MODELISATION A AVEC UTILISATION D'UN CRITERE
#                  DE RE-DECOUPAGE DU PAS DE TEMPS PAR EVENT-DRIVEN
# ......................................................................

MA = LIRE_MAILLAGE(FORMAT="MED", PARTITIONNEUR="PTSCOTCH")


MO = AFFE_MODELE(MAILLAGE=MA, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

ACIER = DEFI_MATERIAU(
    ELAS=_F(E=195000.0, NU=0.3), ECRO_LINE=_F(D_SIGM_EPSI=1930.0, SY=181.0), PRAGER=_F(C=0.0)
)

CM = AFFE_MATERIAU(MAILLAGE=MA, AFFE=_F(TOUT="OUI", MATER=ACIER))

SIGMA_F = DEFI_FONCTION(
    NOM_PARA="INST",
    VALE=(0.0, 0.0, 1.0, 151.2, 2.0, 257.2, 3.0, 0.0),
    PROL_DROITE="EXCLU",
    PROL_GAUCHE="EXCLU",
)

TAU_F = DEFI_FONCTION(
    NOM_PARA="INST",
    VALE=(0.0, 0.0, 1.0, 93.1, 2.0, 33.1, 3.0, 0.0),
    PROL_DROITE="EXCLU",
    PROL_GAUCHE="EXCLU",
)

TRACTION = AFFE_CHAR_MECA(MODELE=MO, FORCE_FACE=_F(GROUP_MA="GAUCHE", FX=-1.0))

CISAIL = AFFE_CHAR_MECA(
    MODELE=MO,
    FORCE_FACE=(
        _F(GROUP_MA="GAUCHE", FY=-1.0),
        _F(GROUP_MA="DROITE", FY=1.0),
        _F(GROUP_MA="HAUT", FX=1.0),
        _F(GROUP_MA="BAS", FX=-1.0),
    ),
)

LIAISON = AFFE_CHAR_CINE(
    MODELE=MO,
    MECA_IMPO=(
        _F(GROUP_NO="NO4", DX=0.0, DY=0.0),
        _F(GROUP_NO="NO8", DX=0.0, DY=0.0, DZ=0.0),
        _F(GROUP_NO="NO2", DX=0.0),
        _F(GROUP_NO="NO6", DX=0.0),
    ),
)

L_INST = DEFI_LIST_REEL(
    DEBUT=0.0,
    INTERVALLE=(
        _F(JUSQU_A=0.1, NOMBRE=1),
        _F(JUSQU_A=0.9, NOMBRE=1),
        _F(JUSQU_A=1.0, NOMBRE=1),
        _F(JUSQU_A=2.0, NOMBRE=1),
        _F(JUSQU_A=3.0, NOMBRE=1),
    ),
)

# gestion manuelle avec event-driven

DEFLIST = DEFI_LIST_INST(
    METHODE="MANUEL",
    DEFI_LIST=_F(LIST_INST=L_INST),
    MODELE=MO,
    ECHEC=_F(
        EVENEMENT="DELTA_GRANDEUR",
        SUBD_NIVEAU=15,
        VALE_REF=0.1e-2,
        NOM_CHAM="VARI_ELGA",
        NOM_CMP="V1",
        GROUP_MA="CUBE",
    ),
)

U = command(
    MODELE=MO,
    CHAM_MATER=CM,
    EXCIT=(
        _F(CHARGE=LIAISON),
        _F(CHARGE=TRACTION, FONC_MULT=SIGMA_F),
        _F(CHARGE=CISAIL, FONC_MULT=TAU_F),
    ),
    COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE"),
    INCREMENT=_F(LIST_INST=DEFLIST),
    NEWTON=_F(MATRICE="TANGENTE", REAC_ITER=3),
    # RECH_LINEAIRE=_F(RHO_MAX=150),
    CONVERGENCE=_F(ITER_GLOB_MAXI=15),
    SOLVEUR=_F(POSTTRAITEMENTS="FORCE"),
)

# gestion manuelle avec critere sur radialite

nbincr = 2

L_INST2 = DEFI_LIST_REEL(
    DEBUT=0.0,
    INTERVALLE=(
        _F(JUSQU_A=0.1, NOMBRE=1),
        _F(JUSQU_A=0.9, NOMBRE=1),
        _F(JUSQU_A=1.0, NOMBRE=1),
        _F(JUSQU_A=2.0, NOMBRE=nbincr),
        _F(JUSQU_A=3.0, NOMBRE=nbincr),
    ),
)

# gestion manuelle avec event-driven

DEFLIST2 = DEFI_LIST_INST(
    METHODE="MANUEL", DEFI_LIST=_F(LIST_INST=L_INST2), ECHEC=_F(SUBD_NIVEAU=25)
)

U2 = command(
    MODELE=MO,
    CHAM_MATER=CM,
    EXCIT=(
        _F(CHARGE=LIAISON),
        _F(CHARGE=TRACTION, FONC_MULT=SIGMA_F),
        _F(CHARGE=CISAIL, FONC_MULT=TAU_F),
    ),
    COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE", RESI_RADI_RELA=0.02),
    INCREMENT=_F(LIST_INST=DEFLIST2),
    NEWTON=_F(MATRICE="TANGENTE", REAC_ITER=3),
    # RECH_LINEAIRE=_F(RHO_MAX=150),
    CONVERGENCE=_F(ITER_GLOB_MAXI=15),
    SOLVEUR=_F(POSTTRAITEMENTS="FORCE"),
)

U3 = command(
    MODELE=MO,
    CHAM_MATER=CM,
    EXCIT=(
        _F(CHARGE=LIAISON),
        _F(CHARGE=TRACTION, FONC_MULT=SIGMA_F),
        _F(CHARGE=CISAIL, FONC_MULT=TAU_F),
    ),
    COMPORTEMENT=_F(RELATION="VMIS_ECMI_LINE", RESI_RADI_RELA=0.02),
    INCREMENT=_F(LIST_INST=DEFLIST2),
    NEWTON=_F(MATRICE="TANGENTE", REAC_ITER=2),
    # RECH_LINEAIRE=_F(RHO_MAX=150),
    CONVERGENCE=_F(ITER_GLOB_MAXI=25),
    SOLVEUR=_F(POSTTRAITEMENTS="FORCE"),
)

U = CALC_CHAMP(
    reuse=U,
    RESULTAT=U,
    CONTRAINTE=("SIGM_ELNO"),
    VARI_INTERNE=("VARI_ELNO"),
    DEFORMATION=("EPSI_ELNO", "EPSP_ELGA"),
)


U2 = CALC_CHAMP(
    reuse=U2,
    RESULTAT=U2,
    CONTRAINTE=("SIGM_ELNO"),
    VARI_INTERNE=("VARI_ELNO"),
    DEFORMATION=("EPSI_ELNO", "EPSP_ELGA"),
)


U3 = CALC_CHAMP(
    reuse=U3,
    RESULTAT=U3,
    CONTRAINTE=("SIGM_ELNO"),
    VARI_INTERNE=("VARI_ELNO"),
    DEFORMATION=("EPSI_ELNO", "EPSP_ELGA"),
)


VARI = CREA_CHAMP(
    TYPE_CHAM="ELNO_VARI_R", OPERATION="EXTR", RESULTAT=U, NOM_CHAM="VARI_ELNO", INST=1.0
)

TEST_RESU(
    CHAM_ELEM=_F(
        GROUP_NO="NO2", NOM_CMP="V1", GROUP_MA="CUBE", CHAM_GD=VARI, VALE_CALC=0.020547265463594
    )
)

TEST_RESU(
    RESU=(
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U,
            NOM_CHAM="SIGM_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="SIXX",
            VALE_CALC=151.200000000,
            VALE_REFE=151.19999999999999,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U,
            NOM_CHAM="SIGM_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="SIXY",
            VALE_CALC=93.100000000,
            VALE_REFE=93.099999999999994,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U,
            NOM_CHAM="EPSI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="EPXX",
            VALE_CALC=0.014829714,
            VALE_REFE=0.014829699999999999,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U,
            NOM_CHAM="EPSI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="EPXY",
            VALE_CALC=0.013601401,
            VALE_REFE=0.0136014,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            POINT=1,
            RESULTAT=U,
            NOM_CHAM="EPSP_ELGA",
            NOM_CMP="EPXX",
            VALE_CALC=0.014054329,
            VALE_REFE=0.014054000000000001,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            POINT=1,
            RESULTAT=U,
            NOM_CHAM="EPSP_ELGA",
            NOM_CMP="EPXY",
            VALE_CALC=0.012980734,
            VALE_REFE=0.012980999999999999,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U,
            NOM_CHAM="VARI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="V1",
            VALE_CALC=0.020547265,
            VALE_REFE=0.020547300000000001,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U,
            NOM_CHAM="VARI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="V2",
            VALE_CALC=1.0,
            VALE_REFE=1.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=2.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U,
            NOM_CHAM="EPSI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="EPXX",
            VALE_CALC=0.035327063,
            VALE_REFE=0.035264999999999998,
            PRECISION=6.0e-3,
            GROUP_MA="CUBE",
        ),
        _F(
            INST=2.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U,
            NOM_CHAM="EPSI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="EPXY",
            VALE_CALC=0.020367031,
            VALE_REFE=0.020471,
            PRECISION=6.0e-3,
            GROUP_MA="CUBE",
        ),
        _F(
            INST=2.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            REFERENCE="AUTRE_ASTER",
            POINT=1,
            RESULTAT=U,
            NOM_CHAM="EPSP_ELGA",
            NOM_CMP="EPXX",
            VALE_CALC=0.034008089,
            VALE_REFE=0.033945999999999997,
            PRECISION=6.0e-3,
            GROUP_MA="CUBE",
        ),
        _F(
            INST=2.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            REFERENCE="AUTRE_ASTER",
            POINT=1,
            RESULTAT=U,
            NOM_CHAM="EPSP_ELGA",
            NOM_CMP="EPXY",
            VALE_CALC=0.020146364,
            VALE_REFE=0.020250000000000001,
            PRECISION=6.0e-3,
            GROUP_MA="CUBE",
        ),
        _F(
            INST=2.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U,
            NOM_CHAM="VARI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="V1",
            VALE_CALC=0.042329286,
            VALE_REFE=0.042328999999999999,
            PRECISION=6.0e-3,
            GROUP_MA="CUBE",
        ),
    )
)

TEST_RESU(
    RESU=(
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U2,
            NOM_CHAM="SIGM_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="SIXX",
            VALE_CALC=151.200000000,
            VALE_REFE=151.19999999999999,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U2,
            NOM_CHAM="SIGM_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="SIXY",
            VALE_CALC=93.100000000,
            VALE_REFE=93.099999999999994,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U2,
            NOM_CHAM="EPSI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="EPXX",
            VALE_CALC=0.014829714,
            VALE_REFE=0.014829699999999999,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U2,
            NOM_CHAM="EPSI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="EPXY",
            VALE_CALC=0.013601401,
            VALE_REFE=0.0136014,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            POINT=1,
            RESULTAT=U2,
            NOM_CHAM="EPSP_ELGA",
            NOM_CMP="EPXX",
            VALE_CALC=0.014054329,
            VALE_REFE=0.014054000000000001,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            POINT=1,
            RESULTAT=U2,
            NOM_CHAM="EPSP_ELGA",
            NOM_CMP="EPXY",
            VALE_CALC=0.012980734,
            VALE_REFE=0.012980999999999999,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U2,
            NOM_CHAM="VARI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="V1",
            VALE_CALC=0.020547265,
            VALE_REFE=0.020547300000000001,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U2,
            NOM_CHAM="VARI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="V2",
            VALE_CALC=1.000000000,
            VALE_REFE=1.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=2.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U2,
            NOM_CHAM="EPSI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="EPXX",
            VALE_CALC=0.035339378,
            VALE_REFE=0.035264999999999998,
            PRECISION=1.0e-2,
            GROUP_MA="CUBE",
        ),
        _F(
            INST=2.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U2,
            NOM_CHAM="EPSI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="EPXY",
            VALE_CALC=0.020321852,
            VALE_REFE=0.020471,
            PRECISION=1.0e-2,
            GROUP_MA="CUBE",
        ),
        _F(
            INST=2.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            REFERENCE="AUTRE_ASTER",
            POINT=1,
            RESULTAT=U2,
            NOM_CHAM="EPSP_ELGA",
            NOM_CMP="EPXX",
            VALE_CALC=0.034020403,
            VALE_REFE=0.033945999999999997,
            PRECISION=1.0e-2,
            GROUP_MA="CUBE",
        ),
        _F(
            INST=2.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            REFERENCE="AUTRE_ASTER",
            POINT=1,
            RESULTAT=U2,
            NOM_CHAM="EPSP_ELGA",
            NOM_CMP="EPXY",
            VALE_CALC=0.020101186,
            VALE_REFE=0.020250000000000001,
            PRECISION=1.0e-2,
            GROUP_MA="CUBE",
        ),
        _F(
            INST=2.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U2,
            NOM_CHAM="VARI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="V1",
            VALE_CALC=0.042329286,
            VALE_REFE=0.042328999999999999,
            PRECISION=1.0e-2,
            GROUP_MA="CUBE",
        ),
    )
)

TEST_RESU(
    RESU=(
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U3,
            NOM_CHAM="SIGM_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="SIXX",
            VALE_CALC=151.200000000,
            VALE_REFE=151.19999999999999,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U3,
            NOM_CHAM="SIGM_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="SIXY",
            VALE_CALC=93.100000000,
            VALE_REFE=93.099999999999994,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U3,
            NOM_CHAM="EPSI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="EPXX",
            VALE_CALC=0.014829714,
            VALE_REFE=0.014829699999999999,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U3,
            NOM_CHAM="EPSI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="EPXY",
            VALE_CALC=0.013601401,
            VALE_REFE=0.0136014,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            POINT=1,
            RESULTAT=U3,
            NOM_CHAM="EPSP_ELGA",
            NOM_CMP="EPXX",
            VALE_CALC=0.014054329,
            VALE_REFE=0.014054000000000001,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            POINT=1,
            RESULTAT=U3,
            NOM_CHAM="EPSP_ELGA",
            NOM_CMP="EPXY",
            VALE_CALC=0.012980734,
            VALE_REFE=0.012980999999999999,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U3,
            NOM_CHAM="VARI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="V1",
            VALE_CALC=0.020547265,
            VALE_REFE=0.020547300000000001,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U3,
            NOM_CHAM="VARI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="V2",
            VALE_CALC=1.000000000,
            VALE_REFE=1.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            GROUP_MA="CUBE",
        ),
        _F(
            INST=2.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U3,
            NOM_CHAM="EPSI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="EPXX",
            VALE_CALC=0.035339378,
            VALE_REFE=0.035264999999999998,
            PRECISION=1.0e-2,
            GROUP_MA="CUBE",
        ),
        _F(
            INST=2.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U3,
            NOM_CHAM="EPSI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="EPXY",
            VALE_CALC=0.020321852,
            VALE_REFE=0.020471,
            PRECISION=1.0e-2,
            GROUP_MA="CUBE",
        ),
        _F(
            INST=2.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            REFERENCE="AUTRE_ASTER",
            POINT=1,
            RESULTAT=U3,
            NOM_CHAM="EPSP_ELGA",
            NOM_CMP="EPXX",
            VALE_CALC=0.034020403,
            VALE_REFE=0.033945999999999997,
            PRECISION=1.0e-2,
            GROUP_MA="CUBE",
        ),
        _F(
            INST=2.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            REFERENCE="AUTRE_ASTER",
            POINT=1,
            RESULTAT=U3,
            NOM_CHAM="EPSP_ELGA",
            NOM_CMP="EPXY",
            VALE_CALC=0.020101186,
            VALE_REFE=0.020250000000000001,
            PRECISION=1.0e-2,
            GROUP_MA="CUBE",
        ),
        _F(
            INST=2.0,
            TOLE_MACHINE=(1.0e-6, 1.0e-9),
            REFERENCE="AUTRE_ASTER",
            RESULTAT=U3,
            NOM_CHAM="VARI_ELNO",
            GROUP_NO="NO2",
            NOM_CMP="V1",
            VALE_CALC=0.042329286,
            VALE_REFE=0.042328999999999999,
            PRECISION=1.0e-2,
            GROUP_MA="CUBE",
        ),
    )
)

FIN()
