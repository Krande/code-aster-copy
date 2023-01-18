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

    liste_instant = [
        champ.getTimeValue(champ.getIndexes()[i]) for i in range(len(champ.getIndexes()))
    ]

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
    EXCIT = args.get("EXCIT")
    # Récupération champs thermiques pour modification du champ mater meca
    AFFE_VARC_TEMP = args.get("AFFE_VARC_TEMP")
    # Calcul thermomecanique
    LIST_INST = args.get("LIST_INST")
    CONVERGENCE = args.get("CONVERGENCE")

    # ------------------------------------------------
    #        PREPARATION OF MATERIALFIELDS
    # ------------------------------------------------

    nb_varc_ther = len(AFFE_VARC_TEMP)

    # Creation of a new MaterialField for each case by copy of intial MaterialField

    material = [None] * nb_varc_ther

    for j in range(nb_varc_ther):
        dfkw = {}
        dfkw["VALE_REF"] = AFFE_VARC_TEMP[j]["VALE_REF"]
        dfkw["EVOL"] = AFFE_VARC_TEMP[j]["EVOL"]
        dfkw["NOM_VARC"] = "TEMP"
        material[j] = MaterialWithAddedExternalStateVariable(CHAM_MATER, dfkw, MODELE.getMesh())
        material[j].build()

    # ------------------------------------------------
    #        CALCULS
    # ------------------------------------------------

    thermeca_dict = NonLinearResultDict("ther_dict")

    for i in range(nb_varc_ther):
        if LIST_INST is None:
            LIST_INST2 = DEFI_LIST_INST(DEFI_LIST=_F(VALE=get_list_inst(AFFE_VARC_TEMP[i]["EVOL"])))
        else:
            LIST_INST2 = LIST_INST

        thermeca_dict[AFFE_VARC_TEMP[i]["NOM_CAS"]] = STAT_NON_LINE(
            MODELE=MODELE,
            CHAM_MATER=material[i],
            CARA_ELEM=CARA_ELEM,
            EXCIT=(_F(CHARGE=EXCIT),),
            INCREMENT=_F(LIST_INST=LIST_INST2),
            CONVERGENCE=CONVERGENCE,
        )
    return thermeca_dict
