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

from ..Messages import UTMESS

from ..Cata.Syntax import _F
from ..CodeCommands import CREA_CHAMP, FORMULE
from .Fracture.raff_xfem_zone import RAFF_XFEM_ZONE


def get_nom_maillage_sdfiss(FISS):
    """retourne le nom du maillage associe au concept FISS"""

    nom_ma = FISS.getMesh()
    return nom_ma


def raff_xfem_ops(self, FISSURE, TYPE, **args):
    """
    Macro RAFF_XFEM
    Calcule un indicateur permettant de caracteriser une zone qui sera raffinee.
    L'indicateur est soit la distance (au fond de fissure pour les fissures,
    a l'interface pour les interfaces), soit un indicateur binaire qui vaut
    1 dans la zone d'interet
    """

    # On importe les definitions des commandes a utiliser dans la macro
    # Le nom de la variable doit etre obligatoirement le nom de la commande

    assert TYPE in ("DISTANCE", "ZONE")

    #  recuperation de la liste des fissures/interfaces
    nbfiss = len(FISSURE)

    # on recupere le concept maillage "associe a la sd"
    MA = FISSURE[0].getMesh()
    nom_ma = get_nom_maillage_sdfiss(FISSURE[0])

    # on verifie que toutes les fissures/interfaces sont rattachees au meme
    # maillage
    for i in range(1, nbfiss):
        nom_ma_i = get_nom_maillage_sdfiss(FISSURE[i])
        if nom_ma_i != nom_ma:
            dict_args = dict(valk=(FISSURE[0].getName(), nom_ma, FISSURE[i].getName(), nom_ma_i))
            UTMESS("F", "XFEM2_10", **dict_args)

    # indicateur de type 'DISTANCE'
    if TYPE == "DISTANCE":
        #  formule distance pour une fissure: -r
        __MDISTF = FORMULE(NOM_PARA=("X1", "X2"), VALE="-1.*sqrt(X1**2+X2**2)")
        #  formule distance pour une interface: -r = -|lsn|
        __MDISTI = FORMULE(NOM_PARA=("X1"), VALE="-1.*sqrt(X1**2)")

        __CERR = [None] * nbfiss
        list_err = []
        list_nom_cmp = []
        for_max = "max("

        for i in range(0, nbfiss):
            fiss = FISSURE[i]

            # recuperation du type de discontinuite :'FISSURE' ou 'INTERFACE'
            # si FISSURE   : l'erreur est la distance au fond de fissure
            # si INTERFACE : l'erreur est la distance a l'interface
            typ_ds = fiss.getDiscontinuityType()

            # extraction des champs level sets
            __CHLN = CREA_CHAMP(
                TYPE_CHAM="NOEU_NEUT_R", OPERATION="EXTR", NOM_CHAM="LNNO", FISSURE=fiss
            )

            if typ_ds == "Crack":
                __CHLTB = CREA_CHAMP(
                    TYPE_CHAM="NOEU_NEUT_R", OPERATION="EXTR", NOM_CHAM="LTNO", FISSURE=fiss
                )

                # on renomme le composante X1 en X2
                __CHLT = CREA_CHAMP(
                    TYPE_CHAM="NOEU_NEUT_R",
                    OPERATION="ASSE",
                    MAILLAGE=MA,
                    ASSE=_F(TOUT="OUI", CHAM_GD=__CHLTB, NOM_CMP="X1", NOM_CMP_RESU="X2"),
                )

            # On affecte à chaque noeud du maillage MA la formule __MDISTF ou
            # __MDISTI
            if typ_ds == "Crack":
                __CHFOR = CREA_CHAMP(
                    TYPE_CHAM="NOEU_NEUT_F",
                    OPERATION="AFFE",
                    MAILLAGE=MA,
                    AFFE=_F(TOUT="OUI", NOM_CMP="X1", VALE_F=__MDISTF),
                )
            elif typ_ds == "Interface":
                __CHFOR = CREA_CHAMP(
                    TYPE_CHAM="NOEU_NEUT_F",
                    OPERATION="AFFE",
                    MAILLAGE=MA,
                    AFFE=_F(TOUT="OUI", NOM_CMP="X1", VALE_F=__MDISTI),
                )

            # on evalue en tout noeud le champ de formules
            if typ_ds == "Crack":
                __CERRB = CREA_CHAMP(
                    TYPE_CHAM="NOEU_NEUT_R",
                    OPERATION="EVAL",
                    CHAM_F=__CHFOR,
                    CHAM_PARA=(__CHLN, __CHLT),
                )

            elif typ_ds == "Interface":
                __CERRB = CREA_CHAMP(
                    TYPE_CHAM="NOEU_NEUT_R", OPERATION="EVAL", CHAM_F=__CHFOR, CHAM_PARA=(__CHLN,)
                )

            # champ d'Erreur de la fissure i
            __CERR[i] = CREA_CHAMP(
                TYPE_CHAM="NOEU_NEUT_R",
                OPERATION="ASSE",
                MAILLAGE=MA,
                ASSE=_F(TOUT="OUI", CHAM_GD=__CERRB, NOM_CMP="X1", NOM_CMP_RESU="X" + str(i + 1)),
            )

            list_err.append(__CERR[i])
            list_nom_cmp.append("X" + str(i + 1))
            for_max = for_max + "X" + str(i + 1) + ","

        # si nbfiss = 1, c'est directement X1
        # si nbfiss > 1 : on prend le max des erreurs de chaque fissure
        for_max = for_max + ")"

        if nbfiss == 1:
            __Erreur = FORMULE(NOM_PARA=(list_nom_cmp), VALE="X1")
        else:
            __Erreur = FORMULE(NOM_PARA=(list_nom_cmp), VALE=for_max)

        # Définition de l'erreur en chaque noeud du maillage
        __CHFORM = CREA_CHAMP(
            TYPE_CHAM="NOEU_NEUT_F",
            OPERATION="AFFE",
            MAILLAGE=MA,
            AFFE=_F(TOUT="OUI", NOM_CMP="X1", VALE_F=__Erreur),
        )

        # champ de sortie
        chamout = CREA_CHAMP(
            TYPE_CHAM="NOEU_NEUT_R", OPERATION="EVAL", CHAM_F=__CHFORM, CHAM_PARA=(list_err)
        )

    # indicateur de type 'ZONE'
    elif TYPE == "ZONE":
        __CERR = [None] * nbfiss
        list_asse = []

        for i in range(0, nbfiss):
            __CERR[i] = RAFF_XFEM_ZONE(FISSURE=FISSURE[i], RAYON=args["RAYON"])
            list_asse.append({"CHAM_GD": __CERR[i], "COEF_R": 1.0, "CUMUL": "OUI", "TOUT": "OUI"})

        # champ de sortie
        chamout = CREA_CHAMP(TYPE_CHAM="CART_NEUT_R", OPERATION="ASSE", MAILLAGE=MA, ASSE=list_asse)

    return chamout
