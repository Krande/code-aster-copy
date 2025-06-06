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

from ..Commands import IMPR_RESU, STAT_NON_LINE, DEFI_LIST_INST
from ..Messages import UTMESS
from ..Cata.Syntax import _F
from ..Objects import ExternalStateVariable, MaterialField, NonLinearResultDict
from .affe_materiau_ops import MaterialWithAddedExternalStateVariable


def get_list_inst(champ):
    """
    Fonction récupérant la liste des instants contenus dans le champ variable de commande
    argument : le champ variable de commande
    sortie : liste contenant les instants
    """

    liste_instant = [champ.getTime(index) for index in champ.getIndexes()]

    return liste_instant


def thermeca_mult_ops(self, **args):
    """
    macro THERMECA_MULT
    """

    args = _F(args)

    # Définition modèle mécanique
    MODELE = args.get("MODELE")
    CARA_ELEM = args.get("CARA_ELEM")
    CHAM_MATER = args.get("CHAM_MATER")
    CHAR_MECA_GLOBAL = args.get("CHAR_MECA_GLOBAL")
    # Récupération champs thermiques pour modification du champ mater meca
    CAS_CHARGE_THER = args.get("CAS_CHARGE_THER")
    # Calcul thermomecanique
    LIST_INST = args.get("LIST_INST")
    CONVERGENCE = args.get("CONVERGENCE")

    # ------------------------------------------------
    #        PREPARATION OF MATERIALFIELDS
    # ------------------------------------------------

    nb_varc_ther = len(CAS_CHARGE_THER)

    # Creation of a new MaterialField for each case by copy of intial MaterialField

    material = [None] * nb_varc_ther
    dict_mater_cas = {}
    dict_evol = {}

    for j in range(nb_varc_ther):
        if CAS_CHARGE_THER[j]["EVOL"].getType() == "EVOL_THER":
            dfkw = {}
            dfkw["VALE_REF"] = CAS_CHARGE_THER[j]["VALE_REF"]
            dfkw["EVOL"] = CAS_CHARGE_THER[j]["EVOL"]
            dfkw["NOM_VARC"] = "TEMP"
            dict_mater_cas[CAS_CHARGE_THER[j]["NOM_CAS"]] = MaterialWithAddedExternalStateVariable(
                CHAM_MATER, dfkw, MODELE.getMesh()
            )
            dict_mater_cas[CAS_CHARGE_THER[j]["NOM_CAS"]].build()
            dict_evol[CAS_CHARGE_THER[j]["NOM_CAS"]] = CAS_CHARGE_THER[j]["EVOL"]

        if CAS_CHARGE_THER[j]["EVOL"].getType() == "EVOL_THER_DICT":
            for nom_cas in CAS_CHARGE_THER[j]["EVOL"].keys():
                dfkw = {}
                dfkw["VALE_REF"] = CAS_CHARGE_THER[j]["VALE_REF"]
                dfkw["EVOL"] = CAS_CHARGE_THER[j]["EVOL"][nom_cas]
                dfkw["NOM_VARC"] = "TEMP"
                dict_mater_cas[nom_cas] = MaterialWithAddedExternalStateVariable(
                    CHAM_MATER, dfkw, MODELE.getMesh()
                )
                dict_mater_cas[nom_cas].build()
                dict_evol[nom_cas] = CAS_CHARGE_THER[j]["EVOL"][nom_cas]

    # ------------------------------------------------
    #        CALCULS
    # ------------------------------------------------

    thermeca_dict = NonLinearResultDict("thermeca_dict")

    for nom_cas in dict_mater_cas.keys():
        if LIST_INST is None:
            LIST_INST2 = DEFI_LIST_INST(DEFI_LIST=_F(VALE=get_list_inst(dict_evol[nom_cas])))
        else:
            LIST_INST2 = LIST_INST

        thermeca_dict[nom_cas] = STAT_NON_LINE(
            MODELE=MODELE,
            CHAM_MATER=dict_mater_cas[nom_cas],
            CARA_ELEM=CARA_ELEM,
            EXCIT=(_F(CHARGE=CHAR_MECA_GLOBAL),),
            INCREMENT=_F(LIST_INST=LIST_INST2),
            CONVERGENCE=CONVERGENCE,
        )
    return thermeca_dict
