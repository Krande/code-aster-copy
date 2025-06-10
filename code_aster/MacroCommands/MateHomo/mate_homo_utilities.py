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

from collections import OrderedDict

from ...Cata.Syntax import _F
from ...Commands import (
    DEFI_CONSTANTE,
    DEFI_MATERIAU,
    AFFE_MATERIAU,
    AFFE_MODELE,
    CREA_RESU,
    CREA_CHAMP,
    DEFI_LIST_REEL,
    AFFE_CHAR_CINE,
    MECA_STATIQUE,
    POST_ELEM,
)
from ...Messages import ASSERT, UTMESS
from ...Objects import Function


def create_empty_dictpara(ls_para):
    tabpara = OrderedDict()
    for para in ls_para:
        tabpara[para] = []
    return tabpara


def get_temp_def_alpha(resu):
    """
    Get the TEMP_DEF_ALPHA parameter from a given SD.
    Return default is 20.0
    """

    temp_def_alpha = 20.0
    mate = resu.getMaterialField()
    for m in mate.getVectorOfMaterial():
        for name in m.getMaterialNames():
            try:
                temp_def_alpha = m.getValueReal(name, "TEMP_DEF_ALPHA")
                break
            except RuntimeError:
                pass

    return temp_def_alpha


def parse_mater_groups(type_homo, ls_affe, varc_name, ls_group_tout):
    """
    Cette fonction sert à :
    * Rajouter une couche de vérification sur les propriétés des matériaux affectées.
    * Rédéfinir ces memes matériaux pour ne conserver que E NU LAMBDA afin de replacer
      les calculs stationnaires par des évolutions. Ceci pour avoir équivalence des résultats.
    * Transformer les affectations TOUT=OUI en affectations sur GROUP_MA=BODY
    """
    # material properties are accessed without "_FO" suffix
    mat1, mat2, mat3 = "ELAS", "THER", "THER_NL"
    mandatory_elas, optional_elas = ["E", "NU"], ["RHO", "ALPHA"]
    mandatory_ther, optional_ther = ["LAMBDA"], ["RHO_CP"]
    need_ther = type_homo in ("MASSIF",)

    affe_mod_mate = []
    affe_mod_calc = []

    f_zero = DEFI_CONSTANTE(VALE=0.0)
    ls_temp_def_alpha = []

    for item in ls_affe:

        mater = item["MATER"]

        if "GROUP_MA" in item:
            new_item_mate = {"GROUP_MA": item["GROUP_MA"]}
            new_item_calc = {"GROUP_MA": item["GROUP_MA"]}
        elif "TOUT" in item:
            new_item_mate = {"GROUP_MA": ls_group_tout}
            new_item_calc = {"GROUP_MA": ls_group_tout}
        else:
            ASSERT(False)

        matNames = mater.getMaterialNames()
        if mat1 not in matNames:
            UTMESS("F", "HOMO1_8", valk=(mat1, mater.getName(), mater.userName))

        if need_ther:
            if mat2 not in matNames and mat3 not in matNames:
                UTMESS("F", "HOMO1_9", valk=(mat2, mat3, mater.getName(), mater.userName))

        keyelas = " ".join([m for m in matNames if m in (mat1,)])
        keyther = " ".join([m for m in matNames if m in (mat2, mat3)])

        f_para = {}
        f_para_temp = {}

        parse_list = {
            keyelas: mandatory_elas + optional_elas,
            keyther: mandatory_ther + optional_ther,
        }

        check_list = mandatory_elas if not need_ther else mandatory_elas + mandatory_ther
        temp_def_alpha_current_mat = None
        for key, lspara in parse_list.items():
            missing_in_at_least_one = []
            for p in lspara:
                func = mater.getFunction(key, p)
                if func:
                    if p in missing_in_at_least_one:
                        UTMESS("F", "HOMO1_12", valk=p)
                    f_para[p] = func
                    try:
                        temp_def_alpha_current_mat = mater.getValueReal(key, "TEMP_DEF_ALPHA")
                        ls_temp_def_alpha.append(temp_def_alpha_current_mat)
                    except RuntimeError:
                        pass
                else:
                    try:
                        v = mater.getValueReal(key, p)
                        UTMESS("F", "HOMO1_13", valk=p)
                    except RuntimeError:
                        pass
                missing_in_at_least_one.append(p)

        for p in check_list:
            if not p in f_para:
                UTMESS("F", "HOMO1_10", valk=(p, mater.getName(), mater.userName))

        for p, fp in f_para.items():
            pro = fp.getProperties()
            if not (pro[0] == "CONSTANT" or (pro[0] == "FONCTION" and pro[2] == varc_name)):
                UTMESS("F", "HOMO1_11", valk=(p, mater.getName(), mater.userName, varc_name))

            if pro[0] == "FONCTION" and pro[2] != "TEMP":
                f_para_temp[p] = Function()
                f_para_temp[p].setParameterName("TEMP")
                f_para_temp[p].setInterpolation(pro[1])
                f_para_temp[p].setExtrapolation(pro[4])
                values = fp.getValues()
                x, y = values[: len(values) // 2], values[len(values) // 2 :]
                f_para_temp[p].setValues(x, y)
            else:
                f_para_temp[p] = fp

        elas_fo_kw = {
            "E": f_para_temp["E"],
            "NU": f_para_temp["NU"],
            "ALPHA": f_zero,
            "TEMP_DEF_ALPHA": 20.0,
        }

        ther_fo_kw = {"LAMBDA": f_para_temp.get("LAMBDA", f_zero), "RHO_CP": f_zero}

        new_item_calc["MATER"] = DEFI_MATERIAU(ELAS_FO=_F(**elas_fo_kw), THER_FO=_F(**ther_fo_kw))

        affe_mod_calc.append(new_item_calc)

        if "ALPHA" in f_para_temp:
            ASSERT(temp_def_alpha_current_mat is not None)
            elas_fo_kw["ALPHA"] = f_para_temp["ALPHA"]
            elas_fo_kw["TEMP_DEF_ALPHA"] = temp_def_alpha_current_mat

        if "RHO" in f_para_temp:
            elas_fo_kw["RHO"] = f_para_temp["RHO"]

        if "RHO_CP" in f_para_temp:
            ther_fo_kw["RHO_CP"] = f_para_temp["RHO_CP"]

        new_item_mate["MATER"] = DEFI_MATERIAU(ELAS_FO=_F(**elas_fo_kw), THER_FO=_F(**ther_fo_kw))

        affe_mod_mate.append(new_item_mate)

    if len(set(ls_temp_def_alpha)) > 1:
        UTMESS("F", "HOMO1_14", valk=("TEMP_DEF_ALPHA",))

    return affe_mod_mate, affe_mod_calc


def prepare_alpha_loads(ls_affe_mod_mate, varc_values):
    """
    Cette fonction sert récuperer les coeff ALPHA affecté sur differents zones
    pour en créer une fonction de INST qui servira dans le affe_char_meca pour
    le calcul des correcteurs de dilatation
    """

    ls_alpha_calc = []

    for item in ls_affe_mod_mate:
        group_ma_affe = item["GROUP_MA"]
        mater = item["MATER"]
        f_alpha_temp = mater.getFunction("ELAS", "ALPHA")
        pro = f_alpha_temp.getProperties()
        f_alpha_time = Function()
        f_alpha_time.setParameterName("INST")
        f_alpha_time.setInterpolation(pro[1])
        f_alpha_time.setExtrapolation(pro[4])
        f_alpha_time.setValues(varc_values, [f_alpha_temp(t) for t in varc_values])
        ls_alpha_calc.append({"GROUP_MA": group_ma_affe, "FONC_ALPHA_TIME": f_alpha_time})

    return ls_alpha_calc


def setup_calcul(type_homo, mesh, ls_group_tout, ls_affe, varc_name, varc_values):

    ls_affe_mod_mate, ls_affe_mod_calc = parse_mater_groups(
        type_homo, ls_affe, varc_name, ls_group_tout
    )

    ls_alpha_calc = prepare_alpha_loads(ls_affe_mod_mate, varc_values)

    MODTH = AFFE_MODELE(
        MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MODELISATION="3D", PHENOMENE="THERMIQUE")
    )

    MODME = AFFE_MODELE(
        MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MODELISATION="3D", PHENOMENE="MECANIQUE")
    )

    EVOLVARC = CREA_RESU(
        OPERATION="AFFE",
        TYPE_RESU="EVOL_VARC",
        NOM_CHAM="TEMP",
        AFFE=[
            _F(
                CHAM_GD=CREA_CHAMP(
                    TYPE_CHAM="NOEU_TEMP_R",
                    OPERATION="AFFE",
                    MAILLAGE=mesh,
                    AFFE=_F(GROUP_MA=ls_group_tout, NOM_CMP="TEMP", VALE=value),
                ),
                INST=value,
            )
            for value in varc_values
        ],
    )

    affevarckw = {
        "GROUP_MA": ls_group_tout,
        "NOM_VARC": "TEMP",
        "EVOL": EVOLVARC,
        "NOM_CHAM": "TEMP",
        "VALE_REF": 20.0,
    }

    CHLOIME = AFFE_MATERIAU(
        MODELE=MODME, AFFE=[_F(**affekw) for affekw in ls_affe_mod_mate], AFFE_VARC=_F(**affevarckw)
    )

    CHMATME = AFFE_MATERIAU(
        MODELE=MODME, AFFE=[_F(**affekw) for affekw in ls_affe_mod_calc], AFFE_VARC=_F(**affevarckw)
    )

    CHMATTH = AFFE_MATERIAU(
        MODELE=MODTH, AFFE=[_F(**affekw) for affekw in ls_affe_mod_calc], AFFE_VARC=_F(**affevarckw)
    )

    L_INST = DEFI_LIST_REEL(VALE=varc_values)

    # Calcul des lois de melange
    # ======================================================================

    LOCK_MECA = AFFE_CHAR_CINE(
        MODELE=MODME, MECA_IMPO=(_F(GROUP_MA=ls_group_tout, DX=0.0, DY=0.0, DZ=0.0))
    )

    DEPLMATE = MECA_STATIQUE(
        MODELE=MODME, CHAM_MATER=CHLOIME, LIST_INST=L_INST, EXCIT=(_F(CHARGE=LOCK_MECA))
    )

    return DEPLMATE, MODME, CHMATME, MODTH, CHMATTH, L_INST, ls_alpha_calc


def cross_work(RESU1, RESU2, INST, ls_group_tout):
    """
    Cette fonction sert à effectuer le calcul du travail croisé des correcteurs au travers
    d'un calcul de l'énergie potentielle d'un correcteur ou de leur combinaison
    """

    ASSERT(RESU1.getMesh() is RESU2.getMesh())
    ASSERT(RESU1.getModel() is RESU2.getModel())
    ASSERT(RESU1.getMaterialField() is RESU2.getMaterialField())
    ASSERT(RESU1.getType() == RESU2.getType())
    epot = 0

    RESU_TYPE = RESU1.getType()
    ASSERT(RESU_TYPE in ("EVOL_ELAS", "EVOL_THER"))
    CH_TYPE = "DEPL" if "ELAS" in RESU_TYPE else "TEMP"

    MOD = RESU1.getModel()
    CHMAT = RESU1.getMaterialField()

    CH1 = CREA_CHAMP(
        RESULTAT=RESU1,
        TYPE_CHAM="NOEU_%s_R" % CH_TYPE,
        OPERATION="EXTR",
        NOM_CHAM=CH_TYPE,
        INST=INST,
    )

    EPOT_CH1 = POST_ELEM(
        CHAM_GD=CH1, INST=INST, MODELE=MOD, CHAM_MATER=CHMAT, ENER_POT=_F(GROUP_MA=ls_group_tout)
    )

    epot_ch1 = abs(sum(EPOT_CH1.EXTR_TABLE().values()["TOTALE"]))

    if RESU1 is RESU2:
        work = 2 * epot_ch1

    else:

        CH2 = CREA_CHAMP(
            RESULTAT=RESU2,
            TYPE_CHAM="NOEU_%s_R" % CH_TYPE,
            OPERATION="EXTR",
            NOM_CHAM=CH_TYPE,
            INST=INST,
        )

        EPOT_CH2 = POST_ELEM(
            CHAM_GD=CH2,
            INST=INST,
            MODELE=MOD,
            CHAM_MATER=CHMAT,
            ENER_POT=_F(GROUP_MA=ls_group_tout),
        )
        epot_ch2 = abs(sum(EPOT_CH2.EXTR_TABLE().values()["TOTALE"]))

        EPOT_SOMME = POST_ELEM(
            CHAM_GD=CH1 + CH2,
            INST=INST,
            MODELE=MOD,
            CHAM_MATER=CHMAT,
            ENER_POT=_F(GROUP_MA=ls_group_tout),
        )
        epot_chsomme = abs(sum(EPOT_SOMME.EXTR_TABLE().values()["TOTALE"]))

        work = epot_chsomme - epot_ch1 - epot_ch2

    return work
