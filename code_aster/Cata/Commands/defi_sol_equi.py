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

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *

DEFI_SOL_EQUI = MACRO(
    nom="DEFI_SOL_EQUI",
    op=OPS("code_aster.MacroCommands.defi_sol_equi_ops.defi_sol_equi_ops"),
    sd_prod=table_sdaster,
    fr=tr("Définition des données de sol pour Miss"),
    reentrant="n",
    regles=(
        ENSEMBLE("TABLE_MATER_ELAS", "TABLE_GEQUI_GMAX", "TABLE_AMOR_EQUI"),
        ENSEMBLE("GROUP_MA_DROITE", "GROUP_MA_GAUCHE"),
    ),
    LIEU_SIGNAL=SIMP(
        statut="f",
        typ="TXM",
        into=("AFFLEURANT", "CHAMP_LIBRE"),
        defaut="AFFLEURANT",
        fr=tr("lieu d'imposition du signal"),
    ),
    CHARGEMENT=SIMP(statut="f", typ="TXM", into=("MONO_APPUI", "ONDE_PLANE"), defaut="MONO_APPUI"),
    b_ONDE=BLOC(
        condition="""equal_to("CHARGEMENT", 'ONDE_PLANE')""",
        regles=(
            UN_PARMI("FONC_SIGNAL", "UNITE_TRAN_INIT", "FONC_SIGNAL_X"),
            ENSEMBLE(
                "FONC_SIGNAL_X",
                "FONC_SIGNAL_Y",
                "FONC_SIGNAL_Z",
                "GROUP_MA_ARRETE_1",
                "GROUP_MA_ARRETE_2",
            ),
            EXCLUS("FONC_SIGNAL_X", "LIAISON"),
        ),
        FONC_SIGNAL=SIMP(
            statut="f", typ=(fonction_sdaster), fr=tr("Signal impose d'accelero horizontal")
        ),
        FONC_SIGNAL_X=SIMP(
            statut="f",
            typ=(fonction_sdaster),
            fr=tr("Signal impose d'accelero dans la direction X"),
        ),
        FONC_SIGNAL_Y=SIMP(
            statut="f",
            typ=(fonction_sdaster),
            fr=tr("Signal impose d'accelero dans la direction Y"),
        ),
        FONC_SIGNAL_Z=SIMP(
            statut="f",
            typ=(fonction_sdaster),
            fr=tr("Signal impose d'accelero dans la direction Z"),
        ),
        TOUT_CHAM=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="OUI"),
        TOUT_ACCE=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="OUI"),
        b_TOUTACCE=BLOC(
            condition="""equal_to("TOUT_ACCE", 'NON')""",
            LIST_COUCHE_ACCE=SIMP(
                statut="f",
                typ=grma,
                validators=NoRepeat(),
                max="**",
                fr=tr("liste des couches pour obtention du signal en accélération"),
            ),
        ),
        # Unite d entree de table de signaux
        UNITE_TRAN_INIT=SIMP(statut="f", typ=UnitType(), inout="in"),
        LIAISON=SIMP(statut="f", typ="TXM", into=("PERIODIQUE", "SANS")),
        MASS_PENA=SIMP(statut="f", typ="R", fr=tr("valeur ponctuelle de masse penalisee")),
        LONG_CARA=SIMP(statut="f", typ="R", fr=tr("valeur de longueur caracteristique")),
        NOM_CMP=SIMP(
            statut="f",
            typ="TXM",
            into=("DX", "DY"),
            defaut="DX",
            fr=tr("sollicitation horizontale"),
        ),
        GROUP_MA_ARRETE_1=SIMP(statut="f", typ=grma),
        GROUP_MA_ARRETE_2=SIMP(statut="f", typ=grma),
        b_MAIL3D=BLOC(
            condition="""exists("FONC_SIGNAL_X")""", MAILLAGE=SIMP(statut="o", typ=maillage_sdaster)
        ),
        b_MAIL2D=BLOC(
            condition="""not exists("FONC_SIGNAL_X")""",
            MAILLAGE=SIMP(statut="f", typ=maillage_sdaster),
        ),
        b_BYRNE=BLOC(
            condition="""not exists("GROUP_MA_ARRETE_1")""",
            CORRECTION=SIMP(
                statut="f",
                typ="TXM",
                into=("BYRNE", "SANS"),
                defaut="SANS",
                fr=tr("utilisation formulation Byrne"),
            ),
            b_CORRECTION=BLOC(
                condition="""equal_to("CORRECTION", 'BYRNE')""",
                COEF_KSI=SIMP(
                    statut="f",
                    typ="R",
                    defaut=0.666666667,
                    fr=tr("facteur Ksi_max par couche pour modele Byrne"),
                ),
            ),
        ),
    ),
    b_MONO=BLOC(
        condition="""equal_to("CHARGEMENT", 'MONO_APPUI')""",
        regles=(
            UN_PARMI("FONC_SIGNAL", "DSP", "FONC_SIGNAL_X"),
            ENSEMBLE(
                "FONC_SIGNAL_X",
                "FONC_SIGNAL_Y",
                "FONC_SIGNAL_Z",
                "GROUP_MA_ARRETE_1",
                "GROUP_MA_ARRETE_2",
            ),
        ),
        FONC_SIGNAL=SIMP(
            statut="f", typ=(fonction_sdaster), fr=tr("Signal impose d'accelero horizontal")
        ),
        FONC_SIGNAL_X=SIMP(
            statut="f",
            typ=(fonction_sdaster),
            fr=tr("Signal impose d'accelero dans la direction X"),
        ),
        FONC_SIGNAL_Y=SIMP(
            statut="f",
            typ=(fonction_sdaster),
            fr=tr("Signal impose d'accelero dans la direction Y"),
        ),
        FONC_SIGNAL_Z=SIMP(
            statut="f",
            typ=(fonction_sdaster),
            fr=tr("Signal impose d'accelero dans la direction Z"),
        ),
        DSP=SIMP(
            statut="f", typ=(fonction_sdaster), fr=tr("DSP Signal impose d'accelero horizontal")
        ),
        b_type_dsp=BLOC(
            condition="""exists('DSP')""",
            DUREE=SIMP(
                statut="o",
                typ="R",
                val_min=0.0,
                fr=tr("durée de la phase forte pour facteur de pic"),
            ),
            UNITE_RESU_DSP=SIMP(statut="o", typ=UnitType(), inout="out"),
        ),
        NOM_CMP=SIMP(
            statut="f",
            typ="TXM",
            into=("DX", "DY"),
            defaut="DX",
            fr=tr("sollicitation horizontale"),
        ),
        TOUT_CHAM=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="NON"),
        TOUT_ACCE=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="OUI"),
        b_TOUTACCE=BLOC(
            condition="""equal_to("TOUT_ACCE", 'NON')""",
            LIST_COUCHE_ACCE=SIMP(
                statut="f",
                typ=grma,
                validators=NoRepeat(),
                max="**",
                fr=tr("liste des couches pour obtention du signal en accélération"),
            ),
        ),
        GROUP_MA_ARRETE_1=SIMP(statut="f", typ=grma),
        GROUP_MA_ARRETE_2=SIMP(statut="f", typ=grma),
        b_MAIL3D=BLOC(
            condition="""exists("FONC_SIGNAL_X")""", MAILLAGE=SIMP(statut="o", typ=maillage_sdaster)
        ),
        b_MAIL2D=BLOC(
            condition="""not exists("FONC_SIGNAL_X")""",
            MAILLAGE=SIMP(statut="f", typ=maillage_sdaster),
        ),
        b_BYRNE=BLOC(
            condition="""not exists("GROUP_MA_ARRETE_1")""",
            CORRECTION=SIMP(
                statut="f",
                typ="TXM",
                into=("BYRNE", "SANS"),
                defaut="SANS",
                fr=tr("utilisation formulation Byrne"),
            ),
            b_CORRECTION=BLOC(
                condition="""equal_to("CORRECTION", 'BYRNE')""",
                COEF_KSI=SIMP(
                    statut="f",
                    typ="R",
                    defaut=0.666666667,
                    fr=tr("facteur Ksi_max par couche pour modele Byrne"),
                ),
            ),
        ),
    ),
    CORR_AMOR=SIMP(
        statut="f",
        typ="TXM",
        into=("OUI", "NON"),
        defaut="NON",
        fr=tr("formulation d'amortissement"),
    ),
    LIST_FREQ_SPEC_OSCI=SIMP(
        statut="f", typ=listr8_sdaster, fr=tr("liste de frequences de spectre d oscillateur")
    ),
    LIST_FREQ=SIMP(statut="f", typ=listr8_sdaster, fr=tr("liste de frequences de calcul")),
    MAILLAGE=SIMP(statut="f", typ=maillage_sdaster),
    GROUP_MA_DROITE=SIMP(statut="o", typ=grma),
    GROUP_MA_GAUCHE=SIMP(statut="o", typ=grma),
    GROUP_MA_SUBSTR=SIMP(statut="o", typ=grma),
    GROUP_MA_COL=SIMP(statut="o", typ=grma),
    COEF_VARI_MATE=SIMP(statut="f", typ="R", defaut=1.0, fr=tr("facteur de variation des modules")),
    COEF_AMPL_ACCE=SIMP(
        statut="f", typ="R", defaut=1.0, fr=tr("facteur sur l'amplitude d'accelero")
    ),
    COEF_GAMMA=SIMP(statut="f", typ="R", defaut=0.65, fr=tr("facteur Gamma_max par couche")),
    NMAX_ITER=SIMP(statut="f", typ="I", defaut=10, val_min=1, fr=tr("nombre d'iterations maximum")),
    RESI_RELA=SIMP(statut="f", typ="R", defaut=0.05, fr=tr("tolerance d'arret des iterations")),
    FREQ_COUP=SIMP(statut="f", typ="R", fr=tr("frequence de coupure de filtrage du signal")),
    FREQ_FILTRE=SIMP(statut="f", typ="R", defaut=0.0, fr=tr("frequence filtrage du signal")),
    INTEGRATION=SIMP(
        statut="f",
        typ="TXM",
        into=("TEMPS", "FREQUENCE"),
        defaut="TEMPS",
        fr=tr("methode integration du signal"),
    ),
    SURF=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="NON"),
    b_SURF=BLOC(
        condition="""equal_to("SURF", 'NON')""",
        # Si SURF='NON' nombre de couches enfoncees
        NIVE_COUCH_ENFO=SIMP(statut="o", typ="I"),
        # Si SURF='NON' nombre de recepteurs decoupant les couches enfoncees
        NB_RECEPTEUR=SIMP(statut="f", typ="I", into=(2, 4), defaut=2),
    ),
    # Unites de sortie
    UNITE_TABLE_RESU=SIMP(statut="f", typ=UnitType(), inout="out"),
    UNITE_RESU_TRAN=SIMP(statut="f", typ=UnitType(), defaut=40, inout="out"),
    UNITE_RESU_SPEC=SIMP(statut="f", typ=UnitType(), defaut=55, inout="out"),
    TABLE_MATER_ELAS=SIMP(statut="f", typ=table_sdaster),
    TABLE_GEQUI_GMAX=SIMP(statut="f", typ=table_sdaster),
    TABLE_AMOR_EQUI=SIMP(statut="f", typ=table_sdaster),
    SEPARATEUR=SIMP(
        statut="f",
        typ="TXM",
        defaut=" ",
        fr=tr("Séparateur des colonnes du tableau (ex : ' ', ';'...)"),
    ),
    LIST_EPSI=SIMP(
        statut="f", typ=listr8_sdaster, fr=tr("liste d'abscisses de distorsion de référence")
    ),
    b_ntabl_mater=BLOC(
        condition="""not exists("TABLE_MATER_ELAS")""",
        MATERIAU=FACT(
            statut="f",
            max="**",
            fr=tr("Définition des matériaux"),
            GAMMA=SIMP(statut="o", typ="R", max="**", fr=tr("Abscisses de distorsion")),
            G_GMAX=SIMP(statut="o", typ="R", max="**", fr=tr("Valeurs de reduction de module G")),
            D=SIMP(statut="o", typ="R", max="**", fr=tr("Valeurs de coefficient d'amortissement")),
        ),
        COUCHE=FACT(
            statut="f",
            max="**",
            fr=tr("Définition des couches"),
            EPAIS=SIMP(statut="o", typ="R", fr=tr("Epaisseur de la couche")),
            GROUP_MA=SIMP(statut="o", typ=grma),
            E=SIMP(statut="o", typ="R", fr=tr("Module d'Young")),
            NU=SIMP(statut="o", typ="R", fr=tr("Coefficient de Poisson")),
            RHO=SIMP(statut="o", typ="R", fr=tr("Masse volumique")),
            AMOR_HYST=SIMP(statut="o", typ="R", fr=tr("Coefficient d'amortissement")),
            NUME_MATE=SIMP(statut="o", typ="I", fr=tr("Numéro du matériau")),
            N1=SIMP(statut="f", typ="R", fr=tr("Valeur essai SPT")),
        ),
    ),
    TITRE=SIMP(statut="f", typ="TXM", fr=tr("Titre de la table produite")),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
)
