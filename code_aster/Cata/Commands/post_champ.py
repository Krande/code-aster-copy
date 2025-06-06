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

# person_in_charge: jacques.pellet at edf.fr

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *


def post_champ_prod(RESULTAT, **args):
    if args.get("__all__"):
        return (resultat_sdaster,)
    if AsType(RESULTAT) is not None:
        return AsType(RESULTAT)
    raise CataError("type de concept resultat non prevu")


# liste des options possibles pour les 4 mots clés EXTR_COQUE, EXTR_TUYAU, EXTR_PMF et MIN_MAX_SP :
liste_option_extr = (
    "EPEQ_ELGA",
    "EPEQ_ELNO",
    "EPSI_ELGA",
    "EPSI_ELNO",
    "SIEF_ELGA",
    "SIEF_ELNO",
    "SIEQ_ELGA",
    "SIEQ_ELNO",
    "SIGM_ELGA",
    "SIGM_ELNO",
    "VARC_ELGA",
    "VARC_ELNO",
    "VARI_ELGA",
    "VARI_ELNO",
    "EPVC_ELGA",
    "EPVC_ELNO",
    "EPME_ELGA",
    "EPME_ELNO",
    "EPSP_ELGA",
    "EPSP_ELNO",
)


POST_CHAMP = OPER(
    nom="POST_CHAMP",
    op=155,
    sd_prod=post_champ_prod,
    reentrant="n",
    fr=tr("extraction de champs sur un sous-point. "),
    regles=(
        UN_PARMI("EXTR_COQUE", "EXTR_TUYAU", "EXTR_PMF", "MIN_MAX_SP", "COQU_EXCENT"),
        EXCLUS(
            "TOUT_ORDRE",
            "NUME_ORDRE",
            "INST",
            "FREQ",
            "NUME_MODE",
            "NOEUD_CMP",
            "LIST_INST",
            "LIST_FREQ",
            "LIST_ORDRE",
            "NOM_CAS",
        ),
        EXCLUS("TOUT", "GROUP_MA"),
    ),
    RESULTAT=SIMP(statut="o", typ=resultat_sdaster, fr=tr("Resultat d'une commande globale")),
    # ====
    # Sélection des numéros d'ordre pour lesquels on fait le calcul :
    # ====
    TOUT_ORDRE=SIMP(statut="f", typ="TXM", into=("OUI",)),
    NUME_ORDRE=SIMP(statut="f", typ="I", validators=NoRepeat(), max="**"),
    NUME_MODE=SIMP(statut="f", typ="I", validators=NoRepeat(), max="**"),
    LIST_ORDRE=SIMP(statut="f", typ=listis_sdaster),
    NOEUD_CMP=SIMP(statut="f", typ="TXM", validators=NoRepeat(), max="**"),
    NOM_CAS=SIMP(statut="f", typ="TXM", validators=NoRepeat(), max="**"),
    INST=SIMP(statut="f", typ="R", validators=NoRepeat(), max="**"),
    LIST_INST=SIMP(statut="f", typ=listr8_sdaster),
    FREQ=SIMP(statut="f", typ="R", validators=NoRepeat(), max="**"),
    LIST_FREQ=SIMP(statut="f", typ=listr8_sdaster),
    b_acce_reel=BLOC(
        condition="""(exists("FREQ"))or(exists("LIST_FREQ"))or(exists("INST"))or(exists("LIST_INST"))""",
        CRITERE=SIMP(statut="f", typ="TXM", defaut="RELATIF", into=("RELATIF", "ABSOLU")),
        b_prec_rela=BLOC(
            condition="""(equal_to("CRITERE", 'RELATIF'))""",
            PRECISION=SIMP(statut="f", typ="R", defaut=1.0e-6),
        ),
        b_prec_abso=BLOC(
            condition="""(equal_to("CRITERE", 'ABSOLU'))""", PRECISION=SIMP(statut="o", typ="R")
        ),
    ),
    # ====
    # Sélection de la zone géométrique:
    # ====
    TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
    GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
    # ====
    # Extraction sur un sous-point d'une coque :
    # ====
    EXTR_COQUE=FACT(
        statut="f",
        max=1,
        fr=tr("extraction sur un sous-point d'une coque"),
        NOM_CHAM=SIMP(
            statut="o", typ="TXM", validators=NoRepeat(), max="**", into=liste_option_extr
        ),
        NUME_COUCHE=SIMP(
            statut="o", typ="I", val_min=1, fr=tr("numero de couche dans l'épaisseur de la coque")
        ),
        NIVE_COUCHE=SIMP(
            statut="o",
            typ="TXM",
            into=("SUP", "INF", "MOY"),
            fr=tr("position dans l'épaisseur de la couche"),
        ),
    ),
    # ====
    # Extraction sur un sous-point d'un tuyau :
    # ====
    EXTR_TUYAU=FACT(
        statut="f",
        max=1,
        fr=tr("extraction sur un sous-point d'un tuyau"),
        NOM_CHAM=SIMP(
            statut="o", typ="TXM", validators=NoRepeat(), max="**", into=liste_option_extr
        ),
        NUME_COUCHE=SIMP(
            statut="o", typ="I", val_min=1, fr=tr("numero de couche dans l'épaisseur du tuyau")
        ),
        NIVE_COUCHE=SIMP(
            statut="o",
            typ="TXM",
            into=("SUP", "INF", "MOY"),
            fr=tr("position dans l'épaisseur de la couche"),
        ),
        ANGLE=SIMP(
            statut="o",
            typ="I",
            val_min=0,
            val_max=360,
            fr=tr("angle de dépouillement pour les tuyaux, en degrés à partir de la génératrice"),
        ),
    ),
    # ====
    # Extraction sur un sous-point d'une poutre multifibre :
    # ====
    EXTR_PMF=FACT(
        statut="f",
        max=1,
        fr=tr("extraction sur un sous-point d'une poutre multifibre"),
        NOM_CHAM=SIMP(
            statut="o", typ="TXM", validators=NoRepeat(), max="**", into=liste_option_extr
        ),
        NUME_FIBRE=SIMP(
            statut="o", typ="I", val_min=1, fr=tr("numéro de la fibre dans la poutre multifibre")
        ),
    ),
    # ====
    # Extraction des min / max sur les sous-points :
    # ====
    MIN_MAX_SP=FACT(
        statut="f",
        max="**",
        fr=tr("extraction du min/max d'une composante pour un champ"),
        regles=(UN_PARMI("NOM_CMP", "NOM_VARI")),
        NOM_CHAM=SIMP(statut="o", typ="TXM", into=liste_option_extr),
        NOM_CMP=SIMP(statut="f", typ="TXM", fr=tr("nom de la composante")),
        NOM_VARI=SIMP(statut="f", typ="TXM", fr=tr("nom de la variable interne")),
        TYPE_MAXI=SIMP(statut="o", typ="TXM", into=("MAXI", "MINI", "MAXI_ABS", "MINI_ABS")),
        NUME_CHAM_RESU=SIMP(
            statut="o",
            typ="I",
            val_min=1,
            val_max=20,
            fr=tr("Numéro du champ produit. Exemple: 6 produit le champ UT06"),
        ),
    ),
    # ====
    # Calcul des efforts des coques "excentrées" sur le feuillet moyen de la coque :
    # ====
    COQU_EXCENT=FACT(
        statut="f",
        max=2,
        fr=tr("Calcul des efforts d'une coque 'excentrée' sur le feuillet moyen de la coque"),
        NOM_CHAM=SIMP(statut="o", typ="TXM", into=("EFGE_ELNO", "EFGE_ELGA")),
        MODI_PLAN=SIMP(statut="o", typ="TXM", into=("OUI",)),
    ),
)
