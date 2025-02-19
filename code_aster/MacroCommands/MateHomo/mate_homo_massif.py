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

import numpy as np

from ...Cata.Syntax import _F
from ...CodeCommands import (
    AFFE_CHAR_CINE,
    AFFE_CHAR_MECA,
    AFFE_CHAR_MECA_F,
    AFFE_CHAR_THER,
    FORMULE,
    MECA_STATIQUE,
    THER_LINEAIRE,
    CALC_CHAMP,
    POST_ELEM,
    CALC_TABLE,
    CREA_TABLE,
)
from ...Messages import ASSERT
from ...Objects import ThermalResultDict, ElasticResultDict

from . import mate_homo_utilities as utilities

# List of all the parameters in the result table

PARAMASSIF = [
    "E_L",
    "E_T",
    "E_N",
    "NU_LT",
    "NU_LN",
    "NU_TN",
    "G_LT",
    "G_LN",
    "G_TN",
    "ALPHA_L",
    "ALPHA_T",
    "ALPHA_N",
    "RHO",
    "LAMBDA_L",
    "LAMBDA_T",
    "LAMBDA_N",
    "RHO_CP",
    "A1111",
    "A2222",
    "A3333",
    "A1122",
    "A1133",
    "A2233",
    "A1212",
    "A2323",
    "A3131",
    "NU_TL",
    "NU_NL",
    "NU_NT",
    "K11",
    "K22",
    "K33",
    "ISOTRANS",
    "TEMP_DEF_ALPHA",
]


def calc_corr_massif_syme(MODME, CHMATME, MODTH, CHMATTH, L_INST, alpha_calc, ls_group_ma):
    """Compute the elastic and thermal correctors for MASSIF (Bulk) case.

    This function performs 7 MECA_STATIQUE and 3 THER_LINEAIRE.
    If the group `face_int` is present performs an additional MECA_STATIQUE to compute
    the internal pressure corrector.

    The computation of the homogeneus parameters for several temperature values is done
    by considering the temperature as a pseudo-time value.

    Arguments
    ---------
        modme (Model): Mechanical model.
        matme (MaterialField): Mechanical material field.
        modth (Model): Thermal model.
        matth (MaterialField): Thermal material field.
        linst (ListOfFloats): List of pseudo-time values (homogeneisation temperature values).
        alpha (list): List of dilatation coefficient as function of pseudo-time (temperature).
        groupma (list[str]): List of groups where properties are prescribed.

    Returns
    -------
        elas (ElasticResultDict): Dict of elastic correctors.
        ther (ThermalResultDict): Dict of thermal correctors.
    """

    # Chargements pour calcul des correcteurs MECANIQUES
    # =======================================================================

    SYME_MECA_XX = AFFE_CHAR_CINE(
        MODELE=MODME,
        MECA_IMPO=(
            _F(GROUP_MA="face_xmin", DX=0.0),
            _F(GROUP_MA="face_ymin", DY=0.0),
            _F(GROUP_MA="face_zmin", DZ=0.0),
            _F(GROUP_MA="face_xmax", DX=0.0),
            _F(GROUP_MA="face_ymax", DY=0.0),
            _F(GROUP_MA="face_zmax", DZ=0.0),
        ),
    )

    ANTI_MECA_12 = AFFE_CHAR_CINE(
        MODELE=MODME,
        MECA_IMPO=(
            _F(GROUP_MA="face_xmin", DY=0.0, DZ=0),
            _F(GROUP_MA="face_ymin", DX=0.0, DZ=0),
            _F(GROUP_MA="face_zmin", DZ=0.0),
            _F(GROUP_MA="face_xmax", DY=0.0, DZ=0),
            _F(GROUP_MA="face_ymax", DX=0.0, DZ=0),
            _F(GROUP_MA="face_zmax", DZ=0.0),
        ),
    )

    ANTI_MECA_31 = AFFE_CHAR_CINE(
        MODELE=MODME,
        MECA_IMPO=(
            _F(GROUP_MA="face_xmin", DZ=0.0, DY=0.0),
            _F(GROUP_MA="face_ymin", DY=0.0),
            _F(GROUP_MA="face_zmin", DX=0.0, DY=0.0),
            _F(GROUP_MA="face_xmax", DZ=0.0, DY=0.0),
            _F(GROUP_MA="face_ymax", DY=0.0),
            _F(GROUP_MA="face_zmax", DX=0.0, DY=0.0),
        ),
    )

    ANTI_MECA_23 = AFFE_CHAR_CINE(
        MODELE=MODME,
        MECA_IMPO=(
            _F(GROUP_MA="face_xmin", DX=0.0),
            _F(GROUP_MA="face_ymin", DZ=0.0, DX=0.0),
            _F(GROUP_MA="face_zmin", DY=0.0, DX=0.0),
            _F(GROUP_MA="face_xmax", DX=0.0),
            _F(GROUP_MA="face_ymax", DZ=0.0, DX=0.0),
            _F(GROUP_MA="face_zmax", DY=0.0, DX=0.0),
        ),
    )

    CHAR11 = AFFE_CHAR_MECA(MODELE=MODME, PRE_EPSI=_F(GROUP_MA=ls_group_ma, EPXX=-1.0))

    CHAR22 = AFFE_CHAR_MECA(MODELE=MODME, PRE_EPSI=_F(GROUP_MA=ls_group_ma, EPYY=-1.0))

    CHAR12 = AFFE_CHAR_MECA(MODELE=MODME, PRE_EPSI=_F(GROUP_MA=ls_group_ma, EPXY=-0.5))

    CHAR33 = AFFE_CHAR_MECA(MODELE=MODME, PRE_EPSI=_F(GROUP_MA=ls_group_ma, EPZZ=-1.0))

    CHAR31 = AFFE_CHAR_MECA(MODELE=MODME, PRE_EPSI=_F(GROUP_MA=ls_group_ma, EPXZ=-0.5))

    CHAR23 = AFFE_CHAR_MECA(MODELE=MODME, PRE_EPSI=_F(GROUP_MA=ls_group_ma, EPYZ=-0.5))

    CHARDIL = AFFE_CHAR_MECA_F(
        MODELE=MODME,
        PRE_EPSI=[
            _F(
                GROUP_MA=item["GROUP_MA"],
                EPXX=FORMULE(
                    VALE="-1.0*ALPHA_DIL(INST)",
                    NOM_PARA=("INST",),
                    ALPHA_DIL=item["FONC_ALPHA_TIME"],
                ),
                EPYY=FORMULE(
                    VALE="-1.0*ALPHA_DIL(INST)",
                    NOM_PARA=("INST",),
                    ALPHA_DIL=item["FONC_ALPHA_TIME"],
                ),
                EPZZ=FORMULE(
                    VALE="-1.0*ALPHA_DIL(INST)",
                    NOM_PARA=("INST",),
                    ALPHA_DIL=item["FONC_ALPHA_TIME"],
                ),
            )
            for item in alpha_calc
        ],
    )

    # Chargements pour calcul des correcteurs THERMIQUES
    # =======================================================================

    SYME_THER_11 = AFFE_CHAR_CINE(
        MODELE=MODTH,
        THER_IMPO=(_F(GROUP_MA="face_xmin", TEMP=0.0), _F(GROUP_MA="face_xmax", TEMP=0.0)),
    )

    SYME_THER_22 = AFFE_CHAR_CINE(
        MODELE=MODTH,
        THER_IMPO=(_F(GROUP_MA="face_ymin", TEMP=0.0), _F(GROUP_MA="face_ymax", TEMP=0.0)),
    )

    SYME_THER_33 = AFFE_CHAR_CINE(
        MODELE=MODTH,
        THER_IMPO=(_F(GROUP_MA="face_zmin", TEMP=0.0), _F(GROUP_MA="face_zmax", TEMP=0.0)),
    )

    CHAR1 = AFFE_CHAR_THER(MODELE=MODTH, PRE_GRAD_TEMP=_F(GROUP_MA=ls_group_ma, FLUX_X=-1.0))

    CHAR2 = AFFE_CHAR_THER(MODELE=MODTH, PRE_GRAD_TEMP=_F(GROUP_MA=ls_group_ma, FLUX_Y=-1.0))

    CHAR3 = AFFE_CHAR_THER(MODELE=MODTH, PRE_GRAD_TEMP=_F(GROUP_MA=ls_group_ma, FLUX_Z=-1.0))

    elas_fields = ElasticResultDict()
    ther_fields = ThermalResultDict()

    # Calcul des correcteurs MECANIQUES de DILATATION
    # ======================================================================

    elas_fields["CORR_DILA"] = MECA_STATIQUE(
        MODELE=MODME,
        CHAM_MATER=CHMATME,
        LIST_INST=L_INST,
        EXCIT=(_F(CHARGE=CHARDIL), _F(CHARGE=SYME_MECA_XX)),
        OPTION="SANS",
    )

    # Calcul des correcteurs MECANIQUES
    # ======================================================================

    elas_fields["CORR_MECA11"] = MECA_STATIQUE(
        MODELE=MODME,
        CHAM_MATER=CHMATME,
        LIST_INST=L_INST,
        EXCIT=(_F(CHARGE=CHAR11), _F(CHARGE=SYME_MECA_XX)),
        OPTION="SANS",
    )

    elas_fields["CORR_MECA22"] = MECA_STATIQUE(
        MODELE=MODME,
        CHAM_MATER=CHMATME,
        LIST_INST=L_INST,
        EXCIT=(_F(CHARGE=CHAR22), _F(CHARGE=SYME_MECA_XX)),
        OPTION="SANS",
    )

    elas_fields["CORR_MECA12"] = MECA_STATIQUE(
        MODELE=MODME,
        CHAM_MATER=CHMATME,
        LIST_INST=L_INST,
        EXCIT=(_F(CHARGE=CHAR12), _F(CHARGE=ANTI_MECA_12)),
        OPTION="SANS",
    )

    elas_fields["CORR_MECA33"] = MECA_STATIQUE(
        MODELE=MODME,
        CHAM_MATER=CHMATME,
        LIST_INST=L_INST,
        EXCIT=(_F(CHARGE=CHAR33), _F(CHARGE=SYME_MECA_XX)),
        OPTION="SANS",
    )

    elas_fields["CORR_MECA31"] = MECA_STATIQUE(
        MODELE=MODME,
        CHAM_MATER=CHMATME,
        LIST_INST=L_INST,
        EXCIT=(_F(CHARGE=CHAR31), _F(CHARGE=ANTI_MECA_31)),
        OPTION="SANS",
    )

    elas_fields["CORR_MECA23"] = MECA_STATIQUE(
        MODELE=MODME,
        CHAM_MATER=CHMATME,
        LIST_INST=L_INST,
        EXCIT=(_F(CHARGE=CHAR23), _F(CHARGE=ANTI_MECA_23)),
        OPTION="SANS",
    )

    # Calcul des correcteurs MECANIQUES de pression interne
    # ======================================================================

    if "face_int" in MODME.getMesh().getGroupsOfCells():
        CHAR_PINT = AFFE_CHAR_MECA(MODELE=MODME, PRES_REP=_F(GROUP_MA="face_int", PRES=1.0))
        elas_fields["CORR_PINT"] = MECA_STATIQUE(
            MODELE=MODME,
            CHAM_MATER=CHMATME,
            LIST_INST=L_INST,
            EXCIT=(_F(CHARGE=CHAR_PINT), _F(CHARGE=SYME_MECA_XX)),
            OPTION="SANS",
        )

    # Calcul des correcteurs THERMIQUES
    # ======================================================================
    ther_fields["CORR_THER11"] = THER_LINEAIRE(
        MODELE=MODTH,
        CHAM_MATER=CHMATTH,
        TYPE_CALCUL="STAT",
        INCREMENT=_F(LIST_INST=L_INST),
        EXCIT=(_F(CHARGE=CHAR1), _F(CHARGE=SYME_THER_11)),
    )

    ther_fields["CORR_THER22"] = THER_LINEAIRE(
        MODELE=MODTH,
        CHAM_MATER=CHMATTH,
        TYPE_CALCUL="STAT",
        INCREMENT=_F(LIST_INST=L_INST),
        EXCIT=(_F(CHARGE=CHAR2), _F(CHARGE=SYME_THER_22)),
    )

    ther_fields["CORR_THER33"] = THER_LINEAIRE(
        MODELE=MODTH,
        CHAM_MATER=CHMATTH,
        TYPE_CALCUL="STAT",
        INCREMENT=_F(LIST_INST=L_INST),
        EXCIT=(_F(CHARGE=CHAR3), _F(CHARGE=SYME_THER_33)),
    )

    return elas_fields, ther_fields


def calc_loimel_massif(DEPLMATE, ls_group_tout):
    """Compute the average value of material parameters on the VER mesh.

    In order to use existing operators (CALC_CHAMP and POST_ELEM) this function works with
    a pseudo-result as input, obtained with a 0-load boundary condition.

    List of computed parameters :
       LAME_1 : first Lamé coefficient
       LAME_2 : second Lamé coefficient
       ALPHA_3K : compression modulus
       RHO : density
       RHO_CP : product of density with specific heat
       LAMBDA_THER : thermal conductivity

    Arguments
    ---------
        deplmate (ElasticResult): Mechanical result from 0-load boundary condition.
        groupma (list[str]): List of groups where properties are prescribed.

    Returns
    -------
        values (dict): average properties values as function of pseudo-time (temperature).
    """

    LAME_1 = FORMULE(NOM_PARA=("E", "NU"), VALE="E*NU/((1+NU)*(1-2*NU))")
    LAME_2 = FORMULE(NOM_PARA=("E", "NU"), VALE="E/(2*(1+NU))")
    ALPHA_3K = FORMULE(NOM_PARA=("E", "NU", "ALPHA"), VALE="ALPHA*E/(1-2*NU)")

    RESUMATE = CALC_CHAMP(
        RESULTAT=DEPLMATE,
        GROUP_MA=ls_group_tout,
        PROPRIETES=("MATE_ELGA",),
        CHAM_UTIL=_F(FORMULE=(LAME_1, LAME_2, ALPHA_3K), NOM_CHAM="MATE_ELGA", NUME_CHAM_RESU=1),
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
            NOM_CMP=("X1", "X2", "X3"),
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
    out["LAME1"], out["LAME2"], out["ALPHA3K"] = [
        LAME_INTE.EXTR_TABLE().values()[key] for key in ("INTE_X1", "INTE_X2", "INTE_X3")
    ]

    out["RHO"], out["RHO_CP"], out["LAMBDA_THER"] = [
        MATE_INTE.EXTR_TABLE().values()[key] for key in ("INTE_RHO", "INTE_RHO_CP", "INTE_LAMBDA")
    ]

    return out


def calc_tabpara_massif(DEPLMATE, volume_ver, ls_group_ma, varc_name, ls_varc, **fields):
    """Compute the homogeneus properties values.

    Arguments
    ---------
        deplmate (ElasticResult): Mechanical result from 0-load boundary condition.
        volumever (float): Volume of VER.
        groupma (list[str]): List of groups where properties are prescribed.
        varcname (str): Name of command variable (TEMP | IRRA).
        varcvalue (list[float]): List of temperature at which parameters are computed.
        **fields (ElasticResultDict, ThermalResultDict): corrector fields.


    Returns
    -------
        A_HOM (list[np.ndarray]): Homogeneus elastic matrix for each temperature value.
        K_HOM (list[np.ndarray]): Homogeneus thermal matrix for each temperature value.
        table (Table): Aster table with all the homonegeus parameters (ready for DEFI_MATERIAU).

    """

    CORR_MECA11 = fields["CORR_MECA11"]
    CORR_MECA22 = fields["CORR_MECA22"]
    CORR_MECA33 = fields["CORR_MECA33"]
    CORR_MECA12 = fields["CORR_MECA12"]
    CORR_MECA31 = fields["CORR_MECA31"]
    CORR_MECA23 = fields["CORR_MECA23"]
    CORR_DILA = fields["CORR_DILA"]
    CORR_THER11 = fields["CORR_THER11"]
    CORR_THER22 = fields["CORR_THER22"]
    CORR_THER33 = fields["CORR_THER33"]

    insts_meca = CORR_MECA11.getAccessParameters()["INST"]
    insts_ther = CORR_THER11.getAccessParameters()["INST"]

    ASSERT(len(insts_meca) == len(insts_ther) == len(ls_varc))

    dictpara = utilities.create_empty_dictpara([varc_name] + PARAMASSIF)
    loimel = calc_loimel_massif(DEPLMATE, ls_group_ma)
    tda = utilities.get_temp_def_alpha(DEPLMATE)
    ls_A_hom = {}
    ls_K_hom = {}

    for i, (inst_meca, inst_ther) in enumerate(zip(insts_meca, insts_ther)):

        ASSERT(inst_meca == inst_ther)

        work_dila_11 = utilities.cross_work(CORR_MECA11, CORR_DILA, inst_meca, ls_group_ma)
        work_dila_22 = utilities.cross_work(CORR_MECA22, CORR_DILA, inst_meca, ls_group_ma)
        work_dila_33 = utilities.cross_work(CORR_MECA33, CORR_DILA, inst_meca, ls_group_ma)

        work_ther_11 = utilities.cross_work(CORR_THER11, CORR_THER11, inst_ther, ls_group_ma)
        work_ther_22 = utilities.cross_work(CORR_THER22, CORR_THER22, inst_ther, ls_group_ma)
        work_ther_33 = utilities.cross_work(CORR_THER33, CORR_THER33, inst_ther, ls_group_ma)

        work_meca_11_11 = utilities.cross_work(CORR_MECA11, CORR_MECA11, inst_meca, ls_group_ma)
        work_meca_22_22 = utilities.cross_work(CORR_MECA22, CORR_MECA22, inst_meca, ls_group_ma)
        work_meca_33_33 = utilities.cross_work(CORR_MECA33, CORR_MECA33, inst_meca, ls_group_ma)

        work_meca_11_22 = utilities.cross_work(CORR_MECA11, CORR_MECA22, inst_meca, ls_group_ma)
        work_meca_11_33 = utilities.cross_work(CORR_MECA11, CORR_MECA33, inst_meca, ls_group_ma)
        work_meca_22_33 = utilities.cross_work(CORR_MECA22, CORR_MECA33, inst_meca, ls_group_ma)

        work_meca_12_12 = utilities.cross_work(CORR_MECA12, CORR_MECA12, inst_meca, ls_group_ma)
        work_meca_23_23 = utilities.cross_work(CORR_MECA23, CORR_MECA23, inst_meca, ls_group_ma)
        work_meca_31_31 = utilities.cross_work(CORR_MECA31, CORR_MECA31, inst_meca, ls_group_ma)

        ########################
        # Matrice homogeneisee
        lambda_ther = loimel["LAMBDA_THER"][i]

        K11_hom = (1 / volume_ver) * (lambda_ther - work_ther_11)
        K22_hom = (1 / volume_ver) * (lambda_ther - work_ther_22)
        K33_hom = (1 / volume_ver) * (lambda_ther - work_ther_33)

        # fmt: off
        K_hom = np.array([[K11_hom, 0,       0      ],
                          [0,       K22_hom, 0      ],
                          [0,       0,       K33_hom]])
        # fmt: on
        ls_K_hom[inst_ther] = K_hom

        lambda_meca = loimel["LAME1"][i]
        mu_meca = loimel["LAME2"][i]

        A1111_hom = (1 / volume_ver) * ((lambda_meca + 2 * mu_meca) - work_meca_11_11)
        A2222_hom = (1 / volume_ver) * ((lambda_meca + 2 * mu_meca) - work_meca_22_22)
        A3333_hom = (1 / volume_ver) * ((lambda_meca + 2 * mu_meca) - work_meca_33_33)

        A1122_hom = (1 / volume_ver) * (lambda_meca - work_meca_11_22)
        A1133_hom = (1 / volume_ver) * (lambda_meca - work_meca_11_33)
        A2233_hom = (1 / volume_ver) * (lambda_meca - work_meca_22_33)

        A1212_hom = (1 / volume_ver) * (mu_meca - work_meca_12_12)
        A2323_hom = (1 / volume_ver) * (mu_meca - work_meca_23_23)
        A3131_hom = (1 / volume_ver) * (mu_meca - work_meca_31_31)

        check_isotrop_trans = abs(round((2 * A1212_hom - A1111_hom + A1122_hom) / A1212_hom, 12))

        # fmt: off
        A_hom = np.array([[A1111_hom, A1122_hom, A1133_hom, 0,         0,         0         ],
                          [A1122_hom, A2222_hom, A2233_hom, 0,         0,         0         ],
                          [A1133_hom, A2233_hom, A3333_hom, 0,         0,         0         ],
                          [0,         0,         0,         A1212_hom, 0,         0         ],
                          [0,         0,         0,         0,         A2323_hom, 0         ],
                          [0,         0,         0,         0,         0,         A3131_hom]])
        # fmt: on
        ls_A_hom[inst_meca] = A_hom

        A_inv = np.linalg.inv(A_hom)
        K_inv = np.linalg.inv(K_hom)

        LAMBDA_L, LAMBDA_T, LAMBDA_N = K_inv.diagonal() ** -1

        E_L, E_T, E_N, G_LT, G_LN, G_TN = A_inv.diagonal() ** -1

        Bdil_11 = (1 / volume_ver) * (loimel["ALPHA3K"][i] - work_dila_11)
        Bdil_22 = (1 / volume_ver) * (loimel["ALPHA3K"][i] - work_dila_22)
        Bdil_33 = (1 / volume_ver) * (loimel["ALPHA3K"][i] - work_dila_33)

        ALPHA_L, ALPHA_T, ALPHA_N = np.dot(A_inv, (Bdil_11, Bdil_22, Bdil_33, 0, 0, 0))[:3]

        NU_TL = -A_inv[0, 1] / A_inv[1, 1]
        NU_NL = -A_inv[0, 2] / A_inv[2, 2]
        NU_NT = -A_inv[1, 2] / A_inv[2, 2]

        NU_LT = -A_inv[1, 0] / A_inv[0, 0]
        NU_LN = -A_inv[2, 0] / A_inv[0, 0]
        NU_TN = -A_inv[2, 1] / A_inv[1, 1]

        RHO = (1 / volume_ver) * loimel["RHO"][i]
        RHO_CP = (1 / volume_ver) * loimel["RHO_CP"][i]

        dictpara["E_L"].append(E_L)
        dictpara["E_T"].append(E_T)
        dictpara["E_N"].append(E_N)
        dictpara["NU_LT"].append(NU_LT)
        dictpara["NU_LN"].append(NU_LN)
        dictpara["NU_TN"].append(NU_TN)
        dictpara["G_LT"].append(G_LT)
        dictpara["G_LN"].append(G_LN)
        dictpara["G_TN"].append(G_TN)
        dictpara["ALPHA_L"].append(ALPHA_L)
        dictpara["ALPHA_T"].append(ALPHA_T)
        dictpara["ALPHA_N"].append(ALPHA_N)
        dictpara["LAMBDA_L"].append(LAMBDA_L)
        dictpara["LAMBDA_T"].append(LAMBDA_T)
        dictpara["LAMBDA_N"].append(LAMBDA_N)
        dictpara["A1111"].append(A1111_hom)
        dictpara["A2222"].append(A2222_hom)
        dictpara["A3333"].append(A3333_hom)
        dictpara["A1122"].append(A1122_hom)
        dictpara["A1133"].append(A1133_hom)
        dictpara["A2233"].append(A2233_hom)
        dictpara["A1212"].append(A1212_hom)
        dictpara["A2323"].append(A2323_hom)
        dictpara["A3131"].append(A3131_hom)
        dictpara["NU_TL"].append(NU_TL)
        dictpara["NU_NL"].append(NU_NL)
        dictpara["NU_NT"].append(NU_NT)
        dictpara["K11"].append(K11_hom)
        dictpara["K22"].append(K22_hom)
        dictpara["K33"].append(K33_hom)
        dictpara["RHO"].append(RHO)
        dictpara["RHO_CP"].append(RHO_CP)
        dictpara["ISOTRANS"].append(check_isotrop_trans)
        dictpara["TEMP_DEF_ALPHA"].append(tda)

    dictpara[varc_name] = ls_varc

    tabpara = CREA_TABLE(LISTE=[_F(PARA=para, LISTE_R=values) for para, values in dictpara.items()])

    return ls_A_hom, ls_K_hom, tabpara
