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

import numpy as np

from ...Cata.Syntax import _F
from ...Commands import (
    AFFE_CHAR_CINE,
    AFFE_CHAR_MECA,
    AFFE_CHAR_MECA_F,
    FORMULE,
    MECA_STATIQUE,
    CALC_CHAMP,
    POST_ELEM,
    CALC_TABLE,
    CREA_TABLE,
)
from ...Messages import ASSERT
from ...Objects import ThermalResultDict, ElasticResultDict

from . import mate_homo_utilities as utilities

# fmt: off
PARAPLAQUE = [
    "MEMB_L", "MEMB_T", "MEMB_LT", "MEMB_G_LT",
    "FLEX_L", "FLEX_T", "FLEX_LT", "FLEX_G_LT",
    "CISA_L", "CISA_T",
    "ALPHA_T", "ALPHA_L",
    "RHO", "TEMP_DEF_ALPHA",
]
# fmt: on


def calc_corr_plaque_syme(MODME, CHMATME, MODTH, CHMATTH, L_INST, ls_group_ma, dir_plaque):

    SYME_MECA_XX_mm = AFFE_CHAR_CINE(
        MODELE=MODME,
        MECA_IMPO=(
            _F(GROUP_MA="face_xmin", DX=0.0),
            _F(GROUP_MA="face_ymin", DY=0.0),
            _F(GROUP_MA="face_zmin", DZ=0.0),
            _F(GROUP_MA="face_xmax", DX=0.0),
            _F(GROUP_MA="face_ymax", DY=0.0),
        ),
    )

    ANTI_MECA_12_mm = AFFE_CHAR_CINE(
        MODELE=MODME,
        MECA_IMPO=(
            _F(GROUP_MA="face_xmin", DY=0.0),
            _F(GROUP_MA="face_ymin", DX=0.0),
            _F(GROUP_MA="face_zmin", DZ=0.0),
            _F(GROUP_MA="face_xmax", DY=0.0),
            _F(GROUP_MA="face_ymax", DX=0.0),
        ),
    )

    SYME_MECA_XX_ff = AFFE_CHAR_CINE(
        MODELE=MODME,
        MECA_IMPO=(
            _F(GROUP_NO="lock_rigi", DZ=0.0),
            _F(GROUP_MA="face_xmin", DX=0.0),
            _F(GROUP_MA="face_ymin", DY=0.0),
            _F(GROUP_MA="face_zmin", DX=0.0, DY=0.0),
            _F(GROUP_MA="face_xmax", DX=0.0),
            _F(GROUP_MA="face_ymax", DY=0.0),
        ),
    )

    ANTI_MECA_12_ff = AFFE_CHAR_CINE(
        MODELE=MODME,
        MECA_IMPO=(
            _F(GROUP_NO="lock_rigi", DZ=0.0),
            _F(GROUP_MA="face_xmin", DY=0.0),
            _F(GROUP_MA="face_ymin", DX=0.0),
            _F(GROUP_MA="face_zmin", DX=0.0, DY=0.0),
            _F(GROUP_MA="face_xmax", DY=0.0),
            _F(GROUP_MA="face_ymax", DX=0.0),
        ),
    )

    LOAD_ff_xx = FORMULE(VALE=f"1.0 * {dir_plaque}", NOM_PARA=[dir_plaque])
    LOAD_ff_xy = FORMULE(VALE=f"0.5 * {dir_plaque}", NOM_PARA=[dir_plaque])

    CHAR11_mm = AFFE_CHAR_MECA(MODELE=MODME, PRE_EPSI=_F(GROUP_MA=ls_group_ma, EPXX=-1.0))

    CHAR22_mm = AFFE_CHAR_MECA(MODELE=MODME, PRE_EPSI=_F(GROUP_MA=ls_group_ma, EPYY=-1.0))

    CHAR12_mm = AFFE_CHAR_MECA(MODELE=MODME, PRE_EPSI=_F(GROUP_MA=ls_group_ma, EPXY=-0.5))

    CHAR11_ff = AFFE_CHAR_MECA_F(MODELE=MODME, PRE_EPSI=_F(GROUP_MA=ls_group_ma, EPXX=LOAD_ff_xx))

    CHAR22_ff = AFFE_CHAR_MECA_F(MODELE=MODME, PRE_EPSI=_F(GROUP_MA=ls_group_ma, EPYY=LOAD_ff_xx))

    CHAR12_ff = AFFE_CHAR_MECA_F(MODELE=MODME, PRE_EPSI=_F(GROUP_MA=ls_group_ma, EPXY=LOAD_ff_xy))

    elas_fields = ElasticResultDict()
    ther_fields = ThermalResultDict()
    # Calcul des correcteurs MECANIQUES
    # ======================================================================

    elas_fields["CORR_MECA11_MEMB"] = MECA_STATIQUE(
        MODELE=MODME,
        CHAM_MATER=CHMATME,
        LIST_INST=L_INST,
        EXCIT=(_F(CHARGE=CHAR11_mm), _F(CHARGE=SYME_MECA_XX_mm)),
    )

    elas_fields["CORR_MECA22_MEMB"] = MECA_STATIQUE(
        MODELE=MODME,
        CHAM_MATER=CHMATME,
        LIST_INST=L_INST,
        EXCIT=(_F(CHARGE=CHAR22_mm), _F(CHARGE=SYME_MECA_XX_mm)),
    )

    elas_fields["CORR_MECA12_MEMB"] = MECA_STATIQUE(
        MODELE=MODME,
        CHAM_MATER=CHMATME,
        LIST_INST=L_INST,
        EXCIT=(_F(CHARGE=CHAR12_mm), _F(CHARGE=ANTI_MECA_12_mm)),
    )

    elas_fields["CORR_MECA11_FLEX"] = MECA_STATIQUE(
        MODELE=MODME,
        CHAM_MATER=CHMATME,
        LIST_INST=L_INST,
        EXCIT=(_F(CHARGE=CHAR11_ff), _F(CHARGE=SYME_MECA_XX_ff)),
    )

    elas_fields["CORR_MECA22_FLEX"] = MECA_STATIQUE(
        MODELE=MODME,
        CHAM_MATER=CHMATME,
        LIST_INST=L_INST,
        EXCIT=(_F(CHARGE=CHAR22_ff), _F(CHARGE=SYME_MECA_XX_ff)),
    )

    elas_fields["CORR_MECA12_FLEX"] = MECA_STATIQUE(
        MODELE=MODME,
        CHAM_MATER=CHMATME,
        LIST_INST=L_INST,
        EXCIT=(_F(CHARGE=CHAR12_ff), _F(CHARGE=ANTI_MECA_12_ff)),
    )

    return elas_fields, ther_fields


def calc_loimel_plaque(DEPLMATE, ls_group_tout, dir_plaque):

    LAME_1_mm = FORMULE(NOM_PARA=("E", "NU"), VALE="E*NU/((1+NU)*(1-2*NU))")
    LAME_2_mm = FORMULE(NOM_PARA=("E", "NU"), VALE="E/(2*(1+NU))")

    LAME_1_ff = FORMULE(
        NOM_PARA=("E", "NU", dir_plaque), VALE="{}**2 * E*NU/((1+NU)*(1-2*NU))".format(dir_plaque)
    )
    LAME_2_ff = FORMULE(
        NOM_PARA=("E", "NU", dir_plaque), VALE="{}**2 * E/(2*(1+NU))".format(dir_plaque)
    )

    RESUMATE = CALC_CHAMP(
        RESULTAT=DEPLMATE,
        GROUP_MA=ls_group_tout,
        PROPRIETES=("MATE_ELGA",),
        CHAM_UTIL=_F(
            FORMULE=(LAME_1_mm, LAME_2_mm, LAME_1_ff, LAME_2_ff),
            NOM_CHAM="MATE_ELGA",
            NUME_CHAM_RESU=1,
        ),
    )

    MATE_INTE = POST_ELEM(
        RESULTAT=RESUMATE,
        MODELE=RESUMATE.getModel(),
        INTEGRALE=_F(
            GROUP_MA=ls_group_tout,
            NOM_CMP=("RHO", "RHO_CP", "LAMBDA"),
            NOM_CHAM="MATE_ELGA",
            TYPE_MAILLE="3D",
        ),
    )

    LAME_INTE = POST_ELEM(
        RESULTAT=RESUMATE,
        MODELE=RESUMATE.getModel(),
        INTEGRALE=_F(
            GROUP_MA=ls_group_tout,
            NOM_CMP=("X1", "X2", "X3", "X4"),
            NOM_CHAM="UT01_ELGA",
            TYPE_MAILLE="3D",
        ),
    )

    if len(ls_group_tout) > 1:
        MATE_INTE = CALC_TABLE(
            reuse=MATE_INTE,
            TABLE=MATE_INTE,
            ACTION=(
                _F(OPERATION="FILTRE", NOM_PARA="GROUP_MA", CRIT_COMP="EQ", VALE_K="UNION_GROUP_MA")
            ),
        )

        LAME_INTE = CALC_TABLE(
            reuse=LAME_INTE,
            TABLE=LAME_INTE,
            ACTION=(
                _F(OPERATION="FILTRE", NOM_PARA="GROUP_MA", CRIT_COMP="EQ", VALE_K="UNION_GROUP_MA")
            ),
        )

    out = {}
    out["LAME1_mm"], out["LAME2_mm"] = [
        LAME_INTE.EXTR_TABLE().values()["INTE_%s" % key] for key in ("X1", "X2")
    ]

    out["LAME1_ff"], out["LAME2_ff"] = [
        LAME_INTE.EXTR_TABLE().values()["INTE_%s" % key] for key in ("X3", "X4")
    ]

    out["RHO"], out["RHO_CP"], out["LAMBDA_THER"] = [
        MATE_INTE.EXTR_TABLE().values()["INTE_%s" % key] for key in ("RHO", "RHO_CP", "LAMBDA")
    ]

    return out


def calc_tabpara_plaque(
    DEPLMATE, volume_ver, ls_group_ma, varc_name, ls_varc, dir_plaque, dirthick, **fields
):

    CORR_MECA11_MEMB = fields["CORR_MECA11_MEMB"]
    CORR_MECA22_MEMB = fields["CORR_MECA22_MEMB"]
    CORR_MECA12_MEMB = fields["CORR_MECA12_MEMB"]
    CORR_MECA11_FLEX = fields["CORR_MECA11_FLEX"]
    CORR_MECA22_FLEX = fields["CORR_MECA22_FLEX"]
    CORR_MECA12_FLEX = fields["CORR_MECA12_FLEX"]

    insts_meca = CORR_MECA11_MEMB.getAccessParameters()["INST"]

    ASSERT(len(insts_meca) == len(ls_varc))

    dictpara = utilities.create_empty_dictpara([varc_name] + PARAPLAQUE)
    loimel = calc_loimel_plaque(DEPLMATE, ls_group_ma, dir_plaque)
    h = dirthick[dir_plaque]
    tda = utilities.get_temp_def_alpha(DEPLMATE)

    # fmt: off
    for i, inst_meca in enumerate(insts_meca):

        work_meca_11_11_mm = utilities.cross_work(
            CORR_MECA11_MEMB, CORR_MECA11_MEMB, inst_meca, ls_group_ma
        )
        work_meca_22_22_mm = utilities.cross_work(
            CORR_MECA22_MEMB, CORR_MECA22_MEMB, inst_meca, ls_group_ma
        )
        work_meca_11_22_mm = utilities.cross_work(
            CORR_MECA11_MEMB, CORR_MECA22_MEMB, inst_meca, ls_group_ma
        )
        work_meca_12_12_mm = utilities.cross_work(
            CORR_MECA12_MEMB, CORR_MECA12_MEMB, inst_meca, ls_group_ma
        )

        work_meca_11_11_ff = utilities.cross_work(
            CORR_MECA11_FLEX, CORR_MECA11_FLEX, inst_meca, ls_group_ma
        )
        work_meca_22_22_ff = utilities.cross_work(
            CORR_MECA22_FLEX, CORR_MECA22_FLEX, inst_meca, ls_group_ma
        )
        work_meca_11_22_ff = utilities.cross_work(
            CORR_MECA11_FLEX, CORR_MECA22_FLEX, inst_meca, ls_group_ma
        )
        work_meca_12_12_ff = utilities.cross_work(
            CORR_MECA12_FLEX, CORR_MECA12_FLEX, inst_meca, ls_group_ma
        )

        lambda_meca_mm = loimel["LAME1_mm"][i]
        mu_meca_mm = loimel["LAME2_mm"][i]

        C1111_hom = (2 * h / volume_ver) * ((lambda_meca_mm + 2 * mu_meca_mm) - work_meca_11_11_mm)
        C2222_hom = (2 * h / volume_ver) * ((lambda_meca_mm + 2 * mu_meca_mm) - work_meca_22_22_mm)
        C1122_hom = (2 * h / volume_ver) * (lambda_meca_mm - work_meca_11_22_mm)
        C1212_hom = (2 * h / volume_ver) * (mu_meca_mm - work_meca_12_12_mm)

        lambda_meca_ff = loimel["LAME1_ff"][i]
        mu_meca_ff = loimel["LAME2_ff"][i]

        D1111_hom = (2 * h / volume_ver) * ((lambda_meca_ff + 2 * mu_meca_ff) - work_meca_11_11_ff)
        D2222_hom = (2 * h / volume_ver) * ((lambda_meca_ff + 2 * mu_meca_ff) - work_meca_22_22_ff)
        D1122_hom = (2 * h / volume_ver) * (lambda_meca_ff - work_meca_11_22_ff)
        D1212_hom = (2 * h / volume_ver) * (mu_meca_ff - work_meca_12_12_ff)

        G11_hom = 0.0
        G22_hom = 0.0

        C_hom = np.array([[C1111_hom, C1122_hom, 0         ],
                          [C1122_hom, C2222_hom, 0         ],
                          [0,         0,         C1212_hom ]])

        D_hom = np.array([[D1111_hom, D1122_hom, 0         ],
                          [D1122_hom, D2222_hom, 0         ],
                          [0,         0,         D1212_hom ]])

        G_hom = np.array([[G11_hom,   0       ],
                          [0,         G22_hom ]])

        RHO = 1 / volume_ver * loimel["RHO"][i]

        ALPHA_T, ALPHA_L = [0, 0]

        dictpara["MEMB_L"].append(C1111_hom)
        dictpara["MEMB_T"].append(C2222_hom)
        dictpara["MEMB_LT"].append(C1122_hom)
        dictpara["MEMB_G_LT"].append(C1212_hom)
        dictpara["FLEX_L"].append(D1111_hom)
        dictpara["FLEX_T"].append(D2222_hom)
        dictpara["FLEX_LT"].append(D1122_hom)
        dictpara["FLEX_G_LT"].append(D1212_hom)
        dictpara["CISA_L"].append(G11_hom)
        dictpara["CISA_T"].append(G22_hom)

        dictpara["ALPHA_L"].append(ALPHA_L)
        dictpara["ALPHA_T"].append(ALPHA_T)
        dictpara["RHO"].append(RHO)
        dictpara["TEMP_DEF_ALPHA"].append(tda)

    # fmt: on
    dictpara[varc_name] = ls_varc

    tabpara = CREA_TABLE(LISTE=[_F(PARA=para, LISTE_R=values) for para, values in dictpara.items()])

    return C_hom, D_hom, G_hom, tabpara
