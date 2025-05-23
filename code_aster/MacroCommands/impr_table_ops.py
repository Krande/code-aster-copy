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

# person_in_charge: mathieu.courtois at edf.fr

import os.path
import re

import libaster
from ..Cata.Syntax import _F
from ..CodeCommands import RECU_FONCTION
from ..Helpers import LogicalUnitFile, ReservedUnitUsed
from ..Utilities.misc import fmtF2PY
from ..Messages import UTMESS


def impr_table_ops(self, FORMAT, TABLE, INFO, **args):
    """
    Macro IMPR_TABLE permettant d'imprimer une table dans un fichier.
    Erreurs<S> dans IMPR_TABLE pour ne pas perdre la base.
    """
    macro = "IMPR_TABLE"

    args = _F(args)

    # On importe les definitions des commandes a utiliser dans la macro
    # Le nom de la variable doit etre obligatoirement le nom de la commande

    # ----------------------------------------------
    # 0. Traitement des arguments, initialisations
    # 0.1. Fichier
    nomfich = LogicalUnitFile.filename_from_unit(args["UNITE"])
    if nomfich and os.path.exists(nomfich) and os.stat(nomfich).st_size != 0:
        if FORMAT == "XMGRACE":
            UTMESS("A", "TABLE0_6", valk=nomfich)

    # 0.2. Création des dictionnaires des FILTRES
    Filtre = []
    if args["FILTRE"]:
        for Fi in args["FILTRE"]:
            dF = Fi.cree_dict_valeurs(Fi.mc_liste)
            for mc in list(dF.keys()):
                if dF[mc] is None:
                    del dF[mc]
            Filtre.append(dF)
    # format pour l'impression des filtres
    form_filtre = "\nFILTRE -> NOM_PARA: %-16s CRIT_COMP: %-4s VALE: %s"

    # 0.3. Création de la liste des tables
    # on conserve la liste même si aujourd'hui, on n'en imprime qu'une à la
    # fois
    ltab = [[TABLE.EXTR_TABLE(), TABLE]]

    # 0.4.1. liste des paramètres à conserver
    nom_para = ltab[0][0].para
    if args["NOM_PARA"]:
        nom_para = args["NOM_PARA"]
    if not type(nom_para) in (list, tuple):
        nom_para = [nom_para]

    # ----------------------------------------------
    # Boucle sur les tables
    for tab, sdtab in ltab:
        # ----- 1. Infos de base
        if INFO == 2:
            print("IMPRESSION DE LA TABLE : %s" % sdtab.getName())

        if args["TITRE"]:
            #    tab.titr = os.linesep.join(args['TITRE'] + (tab.titr, ))
            tab.titr = args["TITRE"] + tab.titr

        # ----- 2. Filtres
        for Fi in Filtre:
            col = getattr(tab, Fi["NOM_PARA"])
            # peu importe le type
            opts = [Fi[k] for k in ("VALE", "VALE_I", "VALE_C", "VALE_K") if k in Fi]
            kargs = {}
            for k in ("CRITERE", "PRECISION"):
                if k in Fi:
                    kargs[k] = Fi[k]
            tab = tab & (getattr(col, Fi["CRIT_COMP"])(*opts, **kargs))
            # trace l'operation dans le titre
            # if FORMAT in ('TABLEAU','ASTER'):
            tab.titr += form_filtre % (
                Fi["NOM_PARA"],
                Fi["CRIT_COMP"],
                " ".join([str(v) for v in opts]),
            )

        # ----- 3. Tris
        if args["TRI"]:
            # une seule occurence de TRI
            T0 = args["TRI"][0]
            dT = T0.cree_dict_valeurs(T0.mc_liste)
            tab.sort(CLES=dT["NOM_PARA"], ORDRE=dT["ORDRE"])

        # ----- 4. Impression
        # vérification des paramètres
        for p in nom_para:
            if not p in tab.para:
                UTMESS("A", "TABLE0_7", valk=p)

        # sélection des paramètres et suppression des colonnes vides
        timp = tab.SansColonneVide(nom_para)

        # passage des mots-clés de mise en forme à la méthode Impr
        kargs = args.copy()
        kargs.update({"FORMAT": FORMAT, "FICHIER": nomfich, "dform": {}})
        # pour l'impression des fonctions
        kfonc = {"FORMAT": FORMAT, "FICHIER": nomfich}

        # 4.1. au format TABLEAU
        if FORMAT == "TABLEAU":
            # surcharge par les formats de l'utilisateur
            kargs["dform"] = {
                "chead": args.get("DEBUT_TABLE"),  # None est traité par Table
                "cfoot": args["FIN_TABLE"],
                "csep": args["SEPARATEUR"],
                "ccom": args["COMMENTAIRE"],
                "ccpara": args["COMM_PARA"],
                "cdeb": args["DEBUT_LIGNE"],
                "cfin": args["FIN_LIGNE"],
            }

        # 4.2. au format AGRAF
        elif FORMAT == "AGRAF":
            kargs["dform"] = {"formR": "%12.5E"}
            kfonc["FORMAT"] = "TABLEAU"

        # 4.3. au format XMGRACE et dérivés
        elif FORMAT == "XMGRACE":
            kargs["dform"] = {"formR": "%.8g"}
            kargs["PILOTE"] = args["PILOTE"]
            kfonc["PILOTE"] = args["PILOTE"]

        # 4.4. format spécifié dans les arguments
        if args["FORMAT_R"]:
            kargs["dform"].update({"formR": fmtF2PY(args["FORMAT_R"])})

        # 4.5. regroupement par paramètre : PAGINATION
        if args["PAGINATION"]:
            l_ppag = args["PAGINATION"]
            if not type(l_ppag) in (list, tuple):
                l_ppag = [l_ppag]
            kargs["PAGINATION"] = [p for p in l_ppag if p in nom_para]
            l_para_err = [p for p in l_ppag if not p in nom_para]
            if len(l_para_err) > 0:
                UTMESS("A", "TABLE0_8", valk=l_para_err)

        with ReservedUnitUsed(args["UNITE"]):
            timp.Impr(**kargs)

        # ----- 5. IMPR_FONCTION='OUI'
        if args["IMPR_FONCTION"] == "OUI":
            # cherche parmi les cellules celles qui contiennent un nom de
            # fonction
            dfon = []
            p_extr = set(["FONCTION", "FONCTION_C"])
            p_extr.intersection_update(timp.para)
            if len(p_extr) > 0:
                # on réduit timp aux colonnes FONCTION et FONCTION_C
                textr = timp.__getitem__(list(p_extr))
                for row in textr:
                    for par, cell in list(row.items()):
                        if type(cell) in (str, str):
                            cell = cell.strip()
                            if libaster.debugJeveuxExists("%-19s.PROL" % cell):
                                dfon.append(["%-19s" % cell, par])
                # impression des fonctions trouvées
                for f, par in dfon:
                    __fonc = RECU_FONCTION(
                        TABLE=sdtab,
                        FILTRE=_F(NOM_PARA=par, VALE_K=f),
                        NOM_PARA_TABL=par,
                        TITRE="Fonction %s" % f,
                    )
                    with ReservedUnitUsed(args["UNITE"]):
                        __fonc.Trace(**kfonc)

    return
