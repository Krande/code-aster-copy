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

DEBUT(CODE="OUI", ERREUR=_F(ALARME="EXCEPTION"))

from zzzz510a_macro import sigelmoy

OPTION = ("SIGM_ELGA", "SIGM_ELMOY")

mamec = LIRE_MAILLAGE()

mamec = MODI_MAILLAGE(reuse=mamec, MAILLAGE=mamec, ORIE_PEAU=_F(GROUP_MA_PEAU="haut"))

eps_plas = (
    0.0,
    0.0002,
    0.0005,
    0.001,
    0.002,
    0.004,
    0.006,
    0.008,
    0.01,
    0.015,
    0.02,
    0.03,
    0.04,
    0.05,
    0.06,
    0.07,
    0.08,
    0.09,
    0.1,
    0.12,
    0.14,
    0.16,
    0.18,
    0.2,
    0.25,
    0.3,
    0.35,
    0.4,
    0.45,
    0.5,
    0.55,
    0.6,
    0.65,
    0.7,
    0.75,
    0.8,
    0.85,
    0.9,
    0.95,
    1.0,
    1.05,
    1.1,
    1.15,
    1.2,
    1.25,
)

sig_mb_20 = (
    577.0,
    577.0,
    577.0,
    577.0,
    577.00,
    577.00,
    577.00,
    577.00,
    594.31,
    623.16,
    646.24,
    686.63,
    729.50,
    752.85,
    773.98,
    789.18,
    800.79,
    812.40,
    822.92,
    841.45,
    857.45,
    871.55,
    884.18,
    895.64,
    920.39,
    941.12,
    959.01,
    974.78,
    988.91,
    1001.72,
    1013.45,
    1024.28,
    1034.35,
    1043.75,
    1052.59,
    1060.92,
    1068.80,
    1076.29,
    1083.43,
    1090.24,
    1096.75,
    1103.00,
    1109.01,
    1114.79,
    1120.36,
)

trac = DEFI_FONCTION(
    NOM_PARA="EPSI",
    ABSCISSE=tuple(
        [sig_mb_20[indi] / 204000.0 + eps_plas[indi] for indi in range(len((eps_plas)))]
    ),
    ORDONNEE=sig_mb_20,
    PROL_DROITE="LINEAIRE",
    PROL_GAUCHE="LINEAIRE",
)

matw = DEFI_MATERIAU(
    ELAS=_F(E=204000.0, NU=0.3),
    TRACTION=_F(SIGM=trac),
    WEIBULL=_F(M=20, VOLU_REFE=1.25e-4, SIGM_REFE=2000, SEUIL_EPSP_CUMU=1e-6),
)

chmat = AFFE_MATERIAU(MAILLAGE=mamec, AFFE=_F(GROUP_MA="s_tout", MATER=matw))

p_int = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 0.0, 1.0, -0.5e4), PROL_DROITE="LINEAIRE")

linst = DEFI_LIST_REEL(DEBUT=0, INTERVALLE=(_F(JUSQU_A=0.17, NOMBRE=10)))

listinst = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=linst), ECHEC=_F(SUBD_NIVEAU=8))

# for modelisation in ("AXIS", "AXIS_SI", "AXIS_INCO_UP", "AXIS_INCO_UPG"):
for modelisation in ("AXIS", "AXIS_SI"):

    momec = AFFE_MODELE(
        MAILLAGE=mamec,
        AFFE=_F(
            GROUP_MA=("s_droite", "s_gauche", "haut"),
            PHENOMENE="MECANIQUE",
            MODELISATION=modelisation,
        ),
    )

    chfmu = CREA_CHAMP(
        OPERATION="AFFE",
        TYPE_CHAM="ELGA_NEUT_F",
        MODELE=momec,
        PROL_ZERO="OUI",
        AFFE=_F(
            GROUP_MA="mgrplasfull",
            NOM_CMP="X1",
            VALE_F=FORMULE(NOM_PARA=("X1", "X2"), VALE="X1*X2"),
        ),
    )

    bloc = AFFE_CHAR_CINE(
        MODELE=momec, MECA_IMPO=(_F(GROUP_MA="gauche", DX=0.0), _F(GROUP_MA="bas", DX=0.0, DY=0.0))
    )

    lunif = AFFE_CHAR_MECA_F(MODELE=momec, LIAISON_UNIF=_F(GROUP_MA="droite", DDL=("DX",)))

    char = AFFE_CHAR_MECA_F(MODELE=momec, PRES_REP=_F(GROUP_MA="haut", PRES=p_int))

    for deformation in ("PETIT", "GDEF_LOG"):

        rmec = STAT_NON_LINE(
            MODELE=momec,
            CHAM_MATER=chmat,
            EXCIT=([{"CHARGE": it} for it in [bloc, char]]),
            COMPORTEMENT=_F(RELATION="VMIS_ISOT_TRAC", DEFORMATION=deformation, GROUP_MA="s_tout"),
            CONVERGENCE=_F(ITER_GLOB_MAXI=40),
            RECH_LINEAIRE=_F(),
            INCREMENT=_F(LIST_INST=listinst),
            SUIVI_DDL=_F(
                NOM_CMP="V1",
                NOM_CHAM="VARI_ELGA",
                GROUP_MA="s_tout",
                EVAL_ELGA="MAX",
                EVAL_CHAM="MAX",
            ),
        )

        sigelmoy(rmec, chfmu, "mgrplasfull")

FIN()
