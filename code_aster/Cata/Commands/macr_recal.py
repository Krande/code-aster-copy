# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

# person_in_charge: thomas.de-soza at edf.fr


from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *


def macr_recal_prod(self, **args):
    if args.get("__all__"):
        return listr8_sdaster

    return listr8_sdaster


MACR_RECAL = MACRO(
    nom="MACR_RECAL",
    op=OPS("code_aster.MacroCommands.macr_recal_ops.macr_recal_ops"),
    sd_prod=macr_recal_prod,
    fr=tr(
        "Réalise le recalage des calculs Aster sur des résultats expérimentaux "
        "ou sur d'autres résultats de calculs"
    ),
    regles=(
        UN_PARMI("PARA_OPTI", "LIST_PARA"),
        PRESENT_PRESENT("PARA_OPTI", "COURBE"),
        PRESENT_PRESENT("LIST_PARA", "RESU_EXP"),
        PRESENT_PRESENT("LIST_PARA", "RESU_CALC"),
        EXCLUS("LIST_POIDS", "COURBE"),
    ),
    UNITE_ESCL=SIMP(statut="o", typ=UnitType(), inout="in"),
    RESU_EXP=SIMP(statut="f", typ=not_checked, max="**"),
    COURBE=FACT(
        statut="f",
        max="**",
        FONC_EXP=SIMP(statut="o", typ=(fonction_sdaster)),
        NOM_FONC_CALC=SIMP(statut="o", typ="TXM"),
        PARA_X=SIMP(statut="o", typ="TXM"),
        PARA_Y=SIMP(statut="o", typ="TXM"),
        POIDS=SIMP(statut="f", typ="R"),
    ),
    RESU_CALC=SIMP(statut="f", typ=not_checked, max="**"),
    LIST_PARA=SIMP(statut="f", typ=not_checked, max="**"),
    PARA_OPTI=FACT(
        statut="f",
        max="**",
        NOM_PARA=SIMP(statut="o", typ="TXM"),
        VALE_INI=SIMP(statut="o", typ="R"),
        VALE_MIN=SIMP(statut="o", typ="R"),
        VALE_MAX=SIMP(statut="o", typ="R"),
    ),
    LIST_POIDS=SIMP(statut="f", typ=not_checked, max="**"),
    UNITE_RESU=SIMP(statut="f", typ=UnitType(), defaut=91, inout="out"),
    PARA_DIFF_FINI=SIMP(statut="f", typ="R", defaut=0.00001),
    GRAPHIQUE=FACT(
        statut="f",
        FORMAT=SIMP(statut="f", typ="TXM", defaut="XMGRACE", into=("XMGRACE", "GNUPLOT")),
        AFFICHAGE=SIMP(
            statut="f",
            typ="TXM",
            defaut="TOUTE_ITERATION",
            into=("TOUTE_ITERATION", "ITERATION_FINALE"),
        ),
        UNITE=SIMP(
            statut="f",
            typ=UnitType(),
            defaut=29,
            inout="out",
            fr=tr("Unité logique définissant le fichier (fort.N) dans lequel on écrit"),
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
        ),
    ),
    # Methodes d'optimisation
    # -----------------------
    METHODE=SIMP(
        statut="f",
        typ="TXM",
        defaut="LEVENBERG",
        into=("LEVENBERG", "FMIN", "FMINBFGS", "FMINNCG", "GENETIQUE", "HYBRIDE"),
    ),
    b_genetique_options=BLOC(
        condition="""equal_to("METHODE", 'GENETIQUE') or equal_to("METHODE", 'HYBRIDE')""",
        NB_PARENTS=SIMP(statut="f", typ="I", defaut=10),
        NB_FILS=SIMP(statut="f", typ="I", defaut=5),
        ECART_TYPE=SIMP(statut="f", typ="R", defaut=1.0),
        GRAINE=SIMP(statut="f", typ="I"),
        ITER_ALGO_GENE=SIMP(statut="f", typ="I", defaut=10),
        RESI_ALGO_GENE=SIMP(statut="f", typ="R", defaut=1.0e-3),
    ),
    # Criteres d'arret globaux
    # -------------------------
    ITER_MAXI=SIMP(
        statut="f", typ="I", defaut=10, fr=tr("Nombre maximum d'iterations d'optimisation")
    ),
    ITER_FONC_MAXI=SIMP(
        statut="f", typ="I", defaut=1000, fr=tr("Nombre maximum d'evaluations de la focntionnelle")
    ),
    RESI_GLOB_RELA=SIMP(
        statut="f", typ="R", defaut=1.0e-3, fr=tr("Critere d'arret sur la valeur du residu")
    ),
    TOLE_PARA=SIMP(
        statut="f", typ="R", defaut=1.0e-8, fr=tr("Critere d'arret sur la valeur des parametres")
    ),
    TOLE_FONC=SIMP(
        statut="f",
        typ="R",
        defaut=1.0e-8,
        fr=tr("Critere d'arret sur la valeur de la fonctionnelle"),
    ),
    # Calculs des gradients
    # ---------------------
    b_gradient=BLOC(
        condition="""equal_to("METHODE", 'FMINBFGS')""",
        GRADIENT=SIMP(
            statut="f",
            typ="TXM",
            defaut="NON_CALCULE",
            into=("NON_CALCULE", "NORMAL", "ADIMENSIONNE"),
        ),
    ),
    b_gradient_fmincg=BLOC(
        condition="""equal_to("METHODE", 'FMINNCG')""",
        GRADIENT=SIMP(statut="o", typ="TXM", into=("NORMAL", "ADIMENSIONNE")),
    ),
    b_gradient_levenberg=BLOC(
        condition="""equal_to("METHODE", 'LEVENBERG')""",
        GRADIENT=SIMP(statut="f", typ="TXM", defaut="NORMAL", into=("NORMAL", "ADIMENSIONNE")),
    ),
    # Mode d'evaluation de l'esclave
    # ------------------------------
    CALCUL_ESCLAVE=FACT(
        statut="d",
        #            regles=(PRESENT_PRESENT('MPI_NBNOEUD','MPI_NBCPU'),),
        LANCEMENT=SIMP(
            statut="f", typ="TXM", defaut="INCLUSION", into=("DISTRIBUTION", "INCLUSION")
        ),
        b_eval_distrib=BLOC(
            condition="""equal_to("LANCEMENT", 'DISTRIBUTION')""",
            UNITE_SUIVI=SIMP(
                statut="f",
                typ=UnitType(),
                defaut=29,
                inout="out",
                fr=tr("Affichage de l'output et/ou error des jobs esclaves dans ce fichier"),
            ),
            MODE=SIMP(statut="f", typ="TXM", into=("INTERACTIF", "BATCH")),
            MEMOIRE=SIMP(
                statut="f", typ="I", fr=tr("Memoire demandee pour les calculs esclaves (Mo)")
            ),
            TEMPS=SIMP(
                statut="f", typ="I", fr=tr("Temps demandee pour les calculs esclaves (secondes)")
            ),
            MPI_NBCPU=SIMP(
                statut="f", typ="I", val_min=1, fr=tr("Nombre de cpu pour les calculs MPI")
            ),
            MPI_NBNOEUD=SIMP(statut="f", typ="I", fr=tr("Nombre de noeuds pour les calculs MPI")),
            CLASSE=SIMP(statut="f", typ="TXM", fr=tr("Classe demandee pour les calculs en batch")),
            NMAX_SIMULT=SIMP(
                statut="f",
                typ="I",
                fr=tr(
                    "Nombre de calculs esclaves lances en parallele en mode distribution (non precise = automatique)"
                ),
            ),
        ),
    ),
    DYNAMIQUE=FACT(
        statut="f",
        MODE_EXP=SIMP(statut="o", typ="TXM"),
        MODE_CALC=SIMP(statut="o", typ="TXM"),
        APPARIEMENT_MANUEL=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
    ),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
)
