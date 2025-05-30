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

from ...Cata.Commons import *
from ...Cata.DataStructure import *
from ...Cata.Syntax import *
from ...Objects import Table
from ...Supervis.ExecuteCommand import ExecuteMacro


def calc_ecrevisse_prod(self, CHARGE_MECA, CHARGE_THER1, CHARGE_THER2, TABLE, DEBIT, **args):
    if args.get("__all__"):
        return ([None], [char_meca], [char_ther], [char_ther], [table_sdaster], [table_sdaster])

    self.type_sdprod(CHARGE_MECA, char_meca)
    self.type_sdprod(CHARGE_THER1, char_ther)
    self.type_sdprod(CHARGE_THER2, char_ther)
    self.type_sdprod(TABLE, table_sdaster)
    self.type_sdprod(DEBIT, table_sdaster)
    return None


CALC_ECREVISSE_CATA = MACRO(
    nom="CALC_ECREVISSE",
    op=OPS("code_aster.MacroCommands.Utils.calc_ecrevisse_ops.calc_ecrevisse_ops"),
    sd_prod=calc_ecrevisse_prod,
    reentrant="n",
    regles=(UN_PARMI("LOGICIEL", "VERSION"),),
    #      CONCEPTS SORTANTS : 2 CHARGEMENTS THERMIQUE + 1 MECANIQUE + 2 TABLES POUR LE POST-TRAITEMENT
    #      ********************************************
    CHARGE_MECA=SIMP(statut="o", typ=CO),
    CHARGE_THER1=SIMP(statut="o", typ=CO),
    CHARGE_THER2=SIMP(statut="o", typ=CO),
    TABLE=SIMP(statut="o", typ=CO),
    DEBIT=SIMP(statut="o", typ=CO),
    #      MODELES MECANIQUES
    #      ********************************************
    MODELE_MECA=SIMP(statut="o", typ=modele_sdaster),
    MODELE_THER=SIMP(statut="o", typ=modele_sdaster),
    #      DONNEES GEOMETRIQUES RELATIVES AUX RESULTATS
    #      ********************************************
    RESULTAT=FACT(
        statut="o",
        min=1,
        max=1,
        MECANIQUE=SIMP(statut="o", typ=resultat_sdaster),
        THERMIQUE=SIMP(statut="o", typ=resultat_sdaster),
        regles=(EXCLUS("NUME_ORDRE", "INST"),),
        NUME_ORDRE=SIMP(statut="f", typ="I", validators=NoRepeat()),
        INST=SIMP(statut="f", typ="R", validators=NoRepeat()),
    ),
    #      DONNEES GEOMETRIQUES RELATIVES A LA FISSURE
    #      *******************************************
    FISSURE=FACT(
        statut="o",
        min=1,
        max="**",
        PREFIXE_FICHIER=SIMP(statut="f", typ="TXM", validators=LongStr(1, 8)),
        GROUP_MA=SIMP(
            statut="o",
            typ=grma,
            validators=NoRepeat(),
            min=2,
            max=2,
            fr=tr("Groupe(s) des noeuds definissant les levres de la fissure"),
        ),
        GROUP_NO_ORIG=SIMP(statut="o", typ=grno, validators=NoRepeat(), min=2, max=2),
        GROUP_NO_EXTR=SIMP(statut="o", typ=grno, validators=NoRepeat(), min=2, max=2),
        ZETA=SIMP(
            statut="o",
            typ="R",
            fr=tr("Coefficient de la perte de charge singuliere a l'entree [zeta]"),
        ),
        RUGOSITE=SIMP(statut="o", typ="R", fr=tr("Rugosite absolu (metres) [eps]")),
        OUVERT_REMANENTE=SIMP(statut="o", typ="R", val_min=0.0, fr=tr("Ouverture remanente")),
        TORTUOSITE=SIMP(
            statut="f",
            typ="R",
            defaut=1.0,
            val_min=0.0,
            val_max=1.0,
            fr=tr("Coefficient de tortuosite de la fissure"),
        ),
        SECTION=SIMP(
            statut="o", typ="TXM", into=("ELLIPSE", "RECTANGLE"), fr=tr("Type de section [is]")
        ),
        b_section_ellipse=BLOC(
            condition="SECTION=='ELLIPSE'",
            fr=tr("Fissure a section elliptique"),
            LISTE_COTES_BL=SIMP(
                statut="f",
                typ="R",
                max="**",
                fr=tr("Liste des cotes des points definissant le petit axe de la section"),
                validators=NoRepeat(),
            ),
            LISTE_VAL_BL=SIMP(
                statut="o",
                typ="R",
                max="**",
                fr=tr("Liste des valeurs des points definissant le petit axe de la section"),
            ),
        ),
        b_section_rectangle=BLOC(
            condition="SECTION=='RECTANGLE'",
            fr=tr("Fissure a section rectangulaire"),
            LISTE_COTES_BL=SIMP(
                statut="f",
                typ="R",
                max="**",
                fr=tr("Liste des cotes des points definissant la largeur de la section"),
                validators=NoRepeat(),
            ),
            LISTE_VAL_BL=SIMP(
                statut="o",
                typ="R",
                max="**",
                fr=tr("Liste des valeurs des points definissant la largeur de la section"),
            ),
        ),
    ),
    #      DONNEES RELATIVES A L"ECOULEMENT
    #      ********************************
    ECOULEMENT=FACT(
        statut="o",
        min=1,
        max=1,
        PRES_ENTREE=SIMP(statut="o", typ="R", fr=tr("Pression de stagnation a l'entree (Pa) [pe]")),
        PRES_SORTIE=SIMP(
            statut="o", typ="R", fr=tr("Pression de stagnation a la sortie (Pa) [ps]")
        ),
        FLUIDE_ENTREE=SIMP(
            statut="o",
            typ="I",
            into=(1, 2, 3, 4, 5, 6),
            fr=tr("Condition du fluide a l'entree [iflow]"),
        ),
        b_condition_1=BLOC(
            condition="FLUIDE_ENTREE==1",
            fr=tr("Eau sous-refroidie ou saturee"),
            TEMP_ENTREE=SIMP(statut="o", typ="R", fr=tr("Temperature a l'entree (degres C) [te]")),
        ),
        b_condition_2=BLOC(
            condition="FLUIDE_ENTREE==2",
            fr=tr("Fluide diphasique"),
            TITR_MASS=SIMP(
                statut="o", typ="R", fr=tr("Titre massique eau vap/eau tot a l'entree [xe]")
            ),
        ),
        b_condition_3=BLOC(
            condition="FLUIDE_ENTREE==3",
            fr=tr("Vapeur saturee ou surchauffee"),
            TEMP_ENTREE=SIMP(statut="o", typ="R", fr=tr("Temperature a l'entree (degres C) [te]")),
        ),
        b_condition_4=BLOC(
            condition="FLUIDE_ENTREE==4",
            fr=tr("Air + vapeur surchauffee"),
            TEMP_ENTREE=SIMP(statut="o", typ="R", fr=tr("Temperature a l'entree (degres C) [te]")),
            PRES_PART=SIMP(
                statut="o", typ="R", fr=tr("Pression partielle air en entree (Pa) [pae]")
            ),
        ),
        b_condition_5=BLOC(
            condition="FLUIDE_ENTREE==5",
            fr=tr("Air + vapeur saturee"),
            TITR_MASS=SIMP(
                statut="o", typ="R", fr=tr("Titre massique eau vap/eau tot a l'entree [xe]")
            ),
            PRES_PART=SIMP(
                statut="o", typ="R", fr=tr("Pression partielle air en entree (Pa) [pae]")
            ),
        ),
        b_condition_6=BLOC(
            condition="FLUIDE_ENTREE==6",
            fr=tr("Air seul"),
            TEMP_ENTREE=SIMP(statut="o", typ="R", fr=tr("Temperature a l'entree (degres C) [te]")),
        ),
    ),
    #      CHOIX DES MODELES
    #      *****************
    MODELE_ECRE=FACT(
        statut="o",
        min=1,
        max=1,
        IVENAC=SIMP(
            statut="f",
            typ="I",
            into=(0, 1),
            defaut=0,
            fr=tr("Calcul ECREVISSE avec prise en compte de la vena contracta"),
        ),
        ECOULEMENT=SIMP(
            statut="o",
            typ="TXM",
            into=("SATURATION", "GELE"),
            fr=tr("Type de modele d'ecoulement diphasique [imod]"),
        ),
        b_ecou_gele=BLOC(
            condition="ECOULEMENT=='GELE'",
            fr=tr("Modele d'ecoulement gele"),
            PRESS_EBULLITION=SIMP(
                statut="o", typ="R", fr=tr("Pression d'ebullition [corrp*psat(t)]")
            ),
        ),
        FROTTEMENT=SIMP(
            statut="o",
            typ="I",
            into=(-4, -3, -2, -1, 0, 1, 2, 3, 4, 11, 12, 13, 14, 21, 22, 23, 24),
            fr=tr("Correlation de frottement [ifrot]"),
        ),
        b_frottement=BLOC(
            condition="FROTTEMENT<0",
            fr=tr("Modele d'ecoulement gele"),
            REYNOLDS_LIM=SIMP(statut="o", typ="R", fr=tr("Coefficient de Reynolds limite [relim]")),
            FROTTEMENT_LIM=SIMP(
                statut="o", typ="R", fr=tr("Coefficient de frottement impose [frtlim]")
            ),
        ),
        TRANSFERT_CHAL=SIMP(
            statut="o",
            typ="I",
            into=(-12, -11, -2, -1, 0, 1, 2, 11, 12),
            fr=tr("Transfert de chaleur [ichal]"),
        ),
        b_transchal=BLOC(
            condition="TRANSFERT_CHAL<0",
            fr=tr("Cas diphasique"),
            XMINCH=SIMP(statut="o", typ="R", fr=tr("Titre massique gazeux min [xminch]")),
            XMAXCH=SIMP(statut="o", typ="R", fr=tr("Titre massique gazeux max [xmaxch]")),
        ),
    ),
    #      DONNEES RELATIVES A LA CONVERGENCE NUMERIQUE
    #      ********************************************
    CONVERGENCE=FACT(
        statut="o",
        min=1,
        max=1,
        KGTEST=SIMP(
            statut="f",
            typ="R",
            val_min=0.0e0,
            val_max=1.0e0,
            defaut=0.5e0,
            fr=tr("Parametre de l'algorithme iteratif [kgtest]"),
        ),
        ITER_GLOB_MAXI=SIMP(
            statut="f",
            typ="I",
            defaut=400,
            fr=tr("Nombre maximum d'iterations de la methode de Newton [itnmax]"),
        ),
        CRIT_CONV_DEBI=SIMP(
            statut="f",
            typ="R",
            val_min=0.0e0,
            val_max=1.0e0,
            defaut=1.0e-5,
            fr=tr("Critere de convergence en debit [precdb]"),
        ),
    ),
    #      GENERAL
    #      *******
    COURBES=SIMP(
        statut="f",
        typ="TXM",
        into=("INTERACTIF", "POSTSCRIPT", "AUCUNE"),
        defaut="AUCUNE",
        fr=tr("Generation eventuelle des courbes"),
    ),
    LOGICIEL=SIMP(statut="f", typ="TXM", validators=LongStr(1, 255)),
    VERSION=SIMP(statut="f", typ="TXM", into=("3.2.2",)),
    ENTETE=SIMP(statut="f", typ="TXM", max="**", defaut="Titre du calcul Ecrevisse"),
    IMPRESSION=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
)


class CalcEcrevisse(ExecuteMacro):
    command_name = "CALC_ECREVISSE"
    command_cata = CALC_ECREVISSE_CATA


CALC_ECREVISSE = CalcEcrevisse.run
