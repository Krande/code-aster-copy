# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *


def calc_srm_prod(self, CHAM_DEFO, **args):
    """Return result type of CALC_SRM."""
    if args.get("__all__"):
        return ([table_sdaster], [None, evol_noli])

    if CHAM_DEFO:
        self.type_sdprod(CHAM_DEFO, evol_noli)

    return table_sdaster


CALC_SRM = MACRO(
    nom="CALC_SRM",
    op=OPS("code_aster.MacroCommands.calc_srm_ops.calc_srm_ops"),
    sd_prod=calc_srm_prod,
    reentrant="n",
    docu="U4.xx.xx",
    fr=tr("Calculer le facteur de sécurité d'une slope par méthode de réduction de la résistance"),
    MODELE=SIMP(statut="o", typ=modele_sdaster),
    CHAM_MATER=SIMP(statut="o", typ=cham_mater),
    # Définition de la zone où s'applique la SRM
    regles=UN_PARMI("TOUT", "GROUP_MA"),
    TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
    GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
    # Configuration du calcul non-linéaire
    EXCIT=FACT(
        statut="f",
        max="**",
        CHARGE=SIMP(statut="o", typ=(char_meca, char_cine_meca)),
        FONC_MULT=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        TYPE_CHARGE=SIMP(
            statut="f", typ="TXM", defaut="FIXE_CSTE", into=("FIXE_CSTE", "SUIV", "DIDI")
        ),
    ),
    INCREMENT=C_INCREMENT("MECANIQUE"),
    CONVERGENCE=C_CONVERGENCE("STAT_NON_LINE"),
    COMPORTEMENT=C_COMPORTEMENT("STAT_NON_LINE"),
    # Configuration de la recherche du FS
    FS=FACT(
        statut="f",
        max=1,
        FS_INIT=SIMP(statut="f", typ="R", defaut=1.0, fr=tr("Facteur de sécurité initial")),
        INCR_INIT=SIMP(statut="f", typ="R", defaut=0.1, fr=tr("Incrément initial de FS")),
        RESI_MAXI=SIMP(
            statut="f", typ="R", defaut=0.01, fr=tr("Résidue du FS après le dernier raffinement")
        ),
        ITER_MAXI=SIMP(statut="f", typ="I", defaut=100, fr=tr("Nombre d'itération maximal")),
        METHODE=SIMP(
            statut="f",
            typ="TXM",
            into=("EXPONENTIELLE", "LINEAIRE"),
            defaut="EXPONENTIELLE",
            fr=tr("Loi de variation de l'incrément du FS"),
        ),
        b_lineaire=BLOC(
            condition="""equal_to("METHODE", "LINEAIRE")""",
            ITER_RAFF_LINE=SIMP(statut="o", typ="I", fr=tr("Nombre de raffinement souhaité")),
        ),
    ),
    # Vérification du résultat et post-traitement
    # VERIF = SIMP(statut = 'f', typ = 'TXM', into = ("NON","OUI"), defaut = "NON",
    #            fr = tr("Vérification du facteur de sécurité par la méthode de BISHOP")),
    CHAM_DEFO=SIMP(
        statut="f",
        typ=CO,
        fr=tr(
            "Champ de déformation plastique cumulée visualisant la surface de glissement lors de la rupture"
        ),
    ),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
)
