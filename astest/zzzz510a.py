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


DEBUT(CODE="OUI")

from zzzz510a_macro import sigelmoy

wb_aster = {
    "group_ma": "CYLINDRE",
    "deformation": "PETIT",
    "option": "SIGM_ELMOY",
    "coef_mult": 2.0,
    "m": 24.0,
    "v_0": 1.25e-4,
    "sigm_cnv": 2600.0,
}

M = LIRE_MAILLAGE(FORMAT="MED")

M = DEFI_GROUP(MAILLAGE=M, CREA_GROUP_NO=(_F(GROUP_MA=("SSUP", "SINF"), NOM="SUNION")))

M = DEFI_GROUP(MAILLAGE=M, DETR_GROUP_NO=_F(NOM=("SSUP", "SINF")))

M = DEFI_GROUP(MAILLAGE=M, CREA_GROUP_NO=(_F(GROUP_MA="SSUP"), _F(GROUP_MA="SINF")))

MO = AFFE_MODELE(AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"), MAILLAGE=M)

SIGY = DEFI_FONCTION(
    PROL_GAUCHE="LINEAIRE",
    NOM_PARA="TEMP",
    VALE=(-150.0, 750.0, -100.0, 700.0, -50.0, 650.0),
    PROL_DROITE="LINEAIRE",
)

SIGREF = DEFI_FONCTION(
    PROL_GAUCHE="LINEAIRE",
    NOM_PARA="TEMP",
    VALE=(-150.0, 2600.0, -50.0, 2800.0),
    PROL_DROITE="LINEAIRE",
)

YOUNG = DEFI_CONSTANTE(VALE=2.0e5)

NUT = DEFI_CONSTANTE(VALE=0.29999999999999999)

ALPHAT = DEFI_CONSTANTE(VALE=0)

ET = DEFI_CONSTANTE(VALE=2000.0)

ACIER = DEFI_MATERIAU(
    ECRO_LINE_FO=_F(D_SIGM_EPSI=ET, SY=SIGY),
    ELAS_FO=_F(NU=NUT, ALPHA=ALPHAT, E=YOUNG, TEMP_DEF_ALPHA=0.0),
    WEIBULL_FO=_F(
        M=wb_aster["m"], VOLU_REFE=wb_aster["v_0"], SIGM_CNV=wb_aster["sigm_cnv"], SIGM_REFE=SIGREF
    ),
)

ZERO = DEFI_CONSTANTE(VALE=0.0)

TRAC50 = DEFI_FONCTION(
    PROL_GAUCHE="EXCLU",
    PROL_DROITE="EXCLU",
    NOM_PARA="INST",
    VALE=(0.0, 0.0, 10.0, 20.35, 20.0, 20.30, 30.0, 20.30, 40.0, 30.525),
)

T_T = DEFI_FONCTION(
    PROL_GAUCHE="CONSTANT",
    NOM_PARA="INST",
    VALE=(0.0, -50.0, 20.0, -50.0, 30.0, -150.0, 40.0, -150.0),
    PROL_DROITE="CONSTANT",
)

TEMP_T = CREA_CHAMP(
    AFFE=_F(NOM_CMP="TEMP", TOUT="OUI", VALE_F=T_T),
    TYPE_CHAM="NOEU_TEMP_F",
    MAILLAGE=M,
    OPERATION="AFFE",
)

L_INST = DEFI_LIST_REEL(INTERVALLE=_F(JUSQU_A=40.0, NOMBRE=4), DEBUT=0.0)

RESU_T = CREA_RESU(
    OPERATION="AFFE",
    TYPE_RESU="EVOL_THER",
    AFFE=(
        _F(NOM_CHAM="TEMP", CHAM_GD=TEMP_T, LIST_INST=L_INST, NUME_FIN=4),
        _F(NOM_CHAM="TEMP", CHAM_GD=TEMP_T, INST=41),
    ),
)

CM = AFFE_MATERIAU(
    AFFE=_F(MATER=ACIER, TOUT="OUI"),
    AFFE_VARC=_F(EVOL=RESU_T, NOM_VARC="TEMP", VALE_REF=0.0, TOUT="OUI"),
    MAILLAGE=M,
)


CHARG = AFFE_CHAR_CINE_F(
    MODELE=MO,
    MECA_IMPO=(
        _F(DY=TRAC50, GROUP_NO="SSUP"),
        _F(DY=ZERO, GROUP_NO="SINF"),
        _F(DX=ZERO, DZ=ZERO, GROUP_NO="P5"),
        _F(DZ=ZERO, GROUP_NO="P6"),
    ),
)

RESU_M = STAT_NON_LINE(
    CHAM_MATER=CM,
    MODELE=MO,
    SOLVEUR=_F(METHODE="MULT_FRONT"),
    COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE"),
    EXCIT=_F(CHARGE=CHARG),
    INCREMENT=_F(NUME_INST_FIN=4, LIST_INST=L_INST),
    NEWTON=_F(MATRICE="TANGENTE", REAC_ITER=1),
)

RESU_M = CALC_CHAMP(reuse=RESU_M, RESULTAT=RESU_M, DEFORMATION="EPSG_ELGA", CONTRAINTE="SIEF_ELNO")

chfmu = CREA_CHAMP(
    OPERATION="AFFE",
    TYPE_CHAM="ELGA_NEUT_F",
    MODELE=MO,
    PROL_ZERO="OUI",
    AFFE=_F(
        GROUP_MA="mgrplasfull", NOM_CMP="X1", VALE_F=FORMULE(NOM_PARA=("X1", "X2"), VALE="X1*X2")
    ),
)

sigelmoy(RESU_M, chfmu, "mgrplasfull")

# TEST union des nodes de group_ma dans un seul group_no
test = CA.TestCase()

no_SSUP = set(M.getNodesFromCells("SSUP"))
no_SINF = set(M.getNodesFromCells("SINF"))
no_UNION = set(M.getNodes("SUNION"))
no_SUP_INF = no_SSUP.union(no_SINF)
test.assertEqual(no_SUP_INF, no_UNION)
test.printSummary()


FIN()
