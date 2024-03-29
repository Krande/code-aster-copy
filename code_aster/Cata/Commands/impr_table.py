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

# person_in_charge: mathieu.courtois at edf.fr

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *

IMPR_TABLE = MACRO(
    nom="IMPR_TABLE",
    op=OPS("code_aster.MacroCommands.impr_table_ops.impr_table_ops"),
    sd_prod=None,
    fr=tr("Impression du contenu d'une table dans un fichier"),
    TABLE=SIMP(statut="o", typ=table_sdaster),
    FORMAT=SIMP(
        statut="f",
        typ="TXM",
        defaut="TABLEAU",
        into=("TABLEAU", "ASTER", "TABLEAU_CROISE", "AGRAF", "XMGRACE", "NUMPY"),
    ),
    b_pilote=BLOC(
        condition="""equal_to("FORMAT", 'XMGRACE')""",
        fr=tr("Mots-clés propres à XMGRACE"),
        PILOTE=SIMP(
            statut="f",
            typ="TXM",
            defaut="",
            into=(
                "",
                "POSTSCRIPT",
                "EPS",
                "MIF",
                "SVG",
                "PNM",
                "PNG",
                "JPEG",
                "PDF",
                "INTERACTIF",
                "INTERACTIF_BG",
            ),
            fr=tr(
                "Pilote de sortie, PNG/JPEG/PDF ne sont pas disponibles sur toutes les installations de xmgrace"
            ),
        ),
        UNITE=SIMP(
            statut="f",
            typ=UnitType(),
            defaut=29,
            inout="out",
            fr=tr("Unité logique définissant le fichier (fort.N) dans lequel on écrit"),
        ),
    ),
    b_unite=BLOC(
        condition="""not equal_to("FORMAT", 'XMGRACE')""",
        UNITE=SIMP(
            statut="f",
            typ=UnitType(),
            defaut=8,
            inout="out",
            fr=tr("Unité logique définissant le fichier (fort.N) dans lequel on écrit"),
        ),
    ),
    FILTRE=FACT(
        statut="f",
        max="**",
        NOM_PARA=SIMP(statut="o", typ="TXM"),
        CRIT_COMP=SIMP(
            statut="f",
            typ="TXM",
            defaut="EQ",
            into=(
                "EQ",
                "LT",
                "GT",
                "NE",
                "LE",
                "GE",
                "VIDE",
                "NON_VIDE",
                "MAXI",
                "MAXI_ABS",
                "MINI",
                "MINI_ABS",
            ),
        ),
        b_vale=BLOC(
            condition="""(is_in("CRIT_COMP", ('EQ','NE','GT','LT','GE','LE')))""",
            regles=(UN_PARMI("VALE", "VALE_I", "VALE_K", "VALE_C"),),
            VALE=SIMP(statut="f", typ="R", max="**"),
            VALE_I=SIMP(statut="f", typ="I", max="**"),
            VALE_C=SIMP(statut="f", typ="C", max="**"),
            VALE_K=SIMP(statut="f", typ="TXM", max="**"),
        ),
        b_crit=BLOC(
            condition="""is_in("CRIT_COMP", ('EQ','NE'))""",
            CRITERE=SIMP(statut="f", typ="TXM", defaut="RELATIF", into=("RELATIF", "ABSOLU")),
            PRECISION=SIMP(statut="f", typ="R", defaut=1.0e-3),
        ),
    ),
    TRI=FACT(
        statut="f",
        NOM_PARA=SIMP(statut="o", typ="TXM", validators=NoRepeat(), max="**"),
        ORDRE=SIMP(
            statut="f",
            typ="TXM",
            defaut="CROISSANT",  # max='**',
            into=("CROISSANT", "DECROISSANT"),
        ),
    ),
    PAGINATION=SIMP(statut="f", typ="TXM", max="**"),
    FORMAT_R=SIMP(statut="f", typ="TXM", defaut="E12.5"),
    NOM_PARA=SIMP(statut="f", typ="TXM", validators=NoRepeat(), max="**"),
    IMPR_FONCTION=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
    # Mise en page du tableau ou du graphique
    b_tableau=BLOC(
        condition="""equal_to("FORMAT", 'TABLEAU')""",
        fr=tr("Mots-clés propres au format Tableau"),
        DEBUT_TABLE=SIMP(
            statut="f",
            typ="TXM",
            fr=tr("Entête avant la table " "(par défaut, une ligne de tirets)"),
        ),
        FIN_TABLE=SIMP(statut="f", typ="TXM", defaut="", fr=tr("Texte après la table")),
        SEPARATEUR=SIMP(
            statut="f",
            typ="TXM",
            defaut=" ",
            fr=tr("Séparateur des colonnes du tableau (ex : ' ', ';'...)"),
        ),
        COMMENTAIRE=SIMP(
            statut="f",
            typ="TXM",
            defaut="#",
            fr=tr("Caractère indiquant au traceur de fonction que la ligne peut etre ignorée"),
        ),
        COMM_PARA=SIMP(
            statut="f",
            typ="TXM",
            defaut="",
            fr=tr("Caractère utilisé pour commentariser la ligne des labels de colonnes"),
        ),
        DEBUT_LIGNE=SIMP(statut="f", typ="TXM", defaut="", fr=tr("Caractère de debut de ligne")),
        FIN_LIGNE=SIMP(statut="f", typ="TXM", defaut="\n", fr=tr("Caractère de fin de ligne")),
    ),
    # mise en forme pour les formats qui passent par Graph
    b_forme=BLOC(
        condition="""equal_to("FORMAT", 'XMGRACE')""",
        fr=tr("Données de mise en forme du graphique"),
        # pour la courbe
        LEGENDE=SIMP(statut="f", typ="TXM", fr=tr("Légende associée à la fonction")),
        STYLE=SIMP(
            statut="f", typ="I", val_min=0, fr=tr("Style de la ligne représentant la fonction")
        ),
        COULEUR=SIMP(statut="f", typ="I", val_min=0, fr=tr("Couleur associée à la fonction")),
        MARQUEUR=SIMP(
            statut="f", typ="I", val_min=0, fr=tr("Type du marqueur associé à la fonction")
        ),
        FREQ_MARQUEUR=SIMP(
            statut="f",
            typ="I",
            defaut=0,
            fr=tr("Fréquence d impression du marqueur associé à la fonction"),
        ),
        # format du graphique
        BORNE_X=SIMP(
            statut="f", typ="R", min=2, max=2, fr=tr("Intervalles de variation des abscisses")
        ),
        BORNE_Y=SIMP(
            statut="f", typ="R", min=2, max=2, fr=tr("Intervalles de variation des ordonnées")
        ),
        ECHELLE_X=SIMP(
            statut="f",
            typ="TXM",
            defaut="LIN",
            into=("LIN", "LOG"),
            fr=tr("Type d'échelle pour les abscisses"),
        ),
        ECHELLE_Y=SIMP(
            statut="f",
            typ="TXM",
            defaut="LIN",
            into=("LIN", "LOG"),
            fr=tr("Type d'échelle pour les ordonnées"),
        ),
        GRILLE_X=SIMP(
            statut="f",
            typ="R",
            max=1,
            val_min=0.0,
            val_min_included=False,
            fr=tr("Pas du quadrillage vertical"),
        ),
        GRILLE_Y=SIMP(
            statut="f",
            typ="R",
            max=1,
            val_min=0.0,
            val_min_included=False,
            fr=tr("Pas du quadrillage horizontal"),
        ),
        LEGENDE_X=SIMP(statut="f", typ="TXM", fr=tr("Légende associée à l'axe des abscisses")),
        LEGENDE_Y=SIMP(statut="f", typ="TXM", fr=tr("Légende associée à l'axe des ordonnées")),
    ),
    TITRE=SIMP(statut="f", typ="TXM"),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
)
