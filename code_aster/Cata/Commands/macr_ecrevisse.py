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

# person_in_charge: marina.bottoni at edf.fr


from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *


from ..Commons.c_comportement import compat_syntax


def macr_ecrevisse_prod(self, TABLE, TEMPER, DEBIT, **args):
    if args.get("__all__"):
        return ([evol_noli], [table_sdaster], [evol_ther], [table_sdaster])

    # On definit ici les concepts produits
    self.type_sdprod(TABLE, table_sdaster)
    self.type_sdprod(TEMPER, evol_ther)
    self.type_sdprod(DEBIT, table_sdaster)
    # concept retourne
    return evol_noli


MACR_ECREVISSE = MACRO(
    nom="MACR_ECREVISSE",
    op=OPS("code_aster.MacroCommands.macr_ecrevisse_ops.macr_ecrevisse_ops"),
    compat_syntax=compat_syntax,
    sd_prod=macr_ecrevisse_prod,
    reentrant="f:ETAT_INIT:EVOL_NOLI",
    fr=tr("Procedure de couplage avec Ecrevisse"),
    reuse=SIMP(statut="c", typ=CO),
    regles=(EXCLUS("TEMPER", "ETAT_INIT"), UN_PARMI("LOGICIEL", "VERSION")),
    #      CONCEPT SORTANT
    #      ********************************************
    TABLE=SIMP(statut="f", typ=CO),
    DEBIT=SIMP(statut="f", typ=CO),
    TEMPER=SIMP(statut="f", typ=CO),
    #      ETAT_INITIAL
    #      ********************************************
    ETAT_INIT=FACT(
        statut="f",
        EVOL_NOLI=SIMP(statut="o", typ=evol_noli),
        EVOL_THER=SIMP(statut="o", typ=evol_ther),
        NUME_ORDRE=SIMP(statut="o", typ="I"),
    ),
    #      MODELES MECANIQUES
    #      ********************************************
    MODELE_MECA=SIMP(statut="o", typ=modele_sdaster),
    MODELE_THER=SIMP(statut="o", typ=modele_sdaster),
    #      DONNEES GEOMETRIQUES RELATIVES A LA FISSURE
    #      *******************************************
    FISSURE=FACT(
        statut="o",
        min=1,
        max="**",
        PREFIXE_FICHIER=SIMP(statut="o", typ="TXM", validators=LongStr(1, 8)),
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
        TORTUOSITE=SIMP(
            statut="f",
            typ="R",
            defaut=1.0,
            val_min=0.0,
            val_max=1.0,
            fr=tr("Coefficient de tortuosite de la fissure"),
        ),
        OUVERT_REMANENTE=SIMP(statut="o", typ="R", val_min=0.0, fr=tr("Ouverture remanente")),
        SECTION=SIMP(
            statut="o", typ="TXM", into=("ELLIPSE", "RECTANGLE"), fr=tr("Type de section [is]")
        ),
        b_section_ellipse=BLOC(
            condition="""equal_to("SECTION", 'ELLIPSE')""",
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
            condition="""equal_to("SECTION", 'RECTANGLE')""",
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
        regles=(
            UN_PARMI("PRES_ENTREE", "PRES_ENTREE_FO"),
            UN_PARMI("PRES_SORTIE", "PRES_SORTIE_FO"),
        ),
        PRES_ENTREE=SIMP(statut="f", typ="R", fr=tr("Pression de stagnation a l'entree (Pa) [pe]")),
        PRES_ENTREE_FO=SIMP(
            statut="f",
            typ=(fonction_sdaster, nappe_sdaster, formule),
            fr=tr("Evolution de la pression de stagnation a l'entree (Pa) [pe]"),
        ),
        PRES_SORTIE=SIMP(
            statut="f", typ="R", fr=tr("Pression de stagnation a la sortie (Pa) [ps]")
        ),
        PRES_SORTIE_FO=SIMP(
            statut="f",
            typ=(fonction_sdaster, nappe_sdaster, formule),
            fr=tr("Evolution de la pression de stagnation a la sortie (Pa) [ps]"),
        ),
        FLUIDE_ENTREE=SIMP(
            statut="o",
            typ="I",
            into=(1, 2, 3, 4, 5, 6),
            fr=tr("Condition du fluide a l'entree [iflow]"),
        ),
        b_condition_1=BLOC(
            condition="""equal_to("FLUIDE_ENTREE", 1)""",
            regles=(UN_PARMI("TEMP_ENTREE", "TEMP_ENTREE_FO")),
            fr=tr("Eau sous-refroidie ou saturee"),
            TEMP_ENTREE=SIMP(statut="f", typ="R", fr=tr("Temperature a l'entree (degres C) [te]")),
            TEMP_ENTREE_FO=SIMP(
                statut="f",
                typ=(fonction_sdaster, nappe_sdaster, formule),
                fr=tr("Evolution de la temperature a l'entree (degres C) [te]"),
            ),
        ),
        b_condition_2=BLOC(
            condition="""equal_to("FLUIDE_ENTREE", 2)""",
            regles=(UN_PARMI("TITR_MASS", "TITR_MASS_FO")),
            fr=tr("Fluide diphasique"),
            TITR_MASS=SIMP(
                statut="f", typ="R", fr=tr("Titre massique eau vap/eau tot a l'entree [xe]")
            ),
            TITR_MASS_FO=SIMP(
                statut="f",
                typ=(fonction_sdaster, nappe_sdaster, formule),
                fr=tr("Evolution du titre massique eau vap/eau tot a l'entree [xe]"),
            ),
        ),
        b_condition_3=BLOC(
            condition="""equal_to("FLUIDE_ENTREE", 3)""",
            regles=(UN_PARMI("TEMP_ENTREE", "TEMP_ENTREE_FO")),
            fr=tr("Vapeur saturee ou surchauffee"),
            TEMP_ENTREE=SIMP(statut="f", typ="R", fr=tr("Temperature a l'entree (degres C) [te]")),
            TEMP_ENTREE_FO=SIMP(
                statut="f",
                typ=(fonction_sdaster, nappe_sdaster, formule),
                fr=tr("Evolution de la temperature a l'entree (degres C) [te]"),
            ),
        ),
        b_condition_4=BLOC(
            condition="""equal_to("FLUIDE_ENTREE", 4)""",
            regles=(
                UN_PARMI("TEMP_ENTREE", "TEMP_ENTREE_FO"),
                UN_PARMI("PRES_PART", "PRES_PART_FO"),
            ),
            fr=tr("Air + vapeur surchauffee"),
            TEMP_ENTREE=SIMP(statut="f", typ="R", fr=tr("Temperature a l'entree (degres C) [te]")),
            TEMP_ENTREE_FO=SIMP(
                statut="f",
                typ=(fonction_sdaster, nappe_sdaster, formule),
                fr=tr("Evolution de la temperature a l'entree (degres C) [te]"),
            ),
            PRES_PART=SIMP(
                statut="f", typ="R", fr=tr("Pression partielle air en entree (Pa) [pae]")
            ),
            PRES_PART_FO=SIMP(
                statut="f",
                typ=(fonction_sdaster, nappe_sdaster, formule),
                fr=tr("Evolution de la pression partielle air en entree (Pa) [pae]"),
            ),
        ),
        b_condition_5=BLOC(
            condition="""equal_to("FLUIDE_ENTREE", 5)""",
            regles=(UN_PARMI("TITR_MASS", "TITR_MASS_FO"), UN_PARMI("PRES_PART", "PRES_PART_FO")),
            fr=tr("Air + vapeur saturee"),
            TITR_MASS=SIMP(
                statut="f", typ="R", fr=tr("Titre massique eau vap/eau tot a l'entree [xe]")
            ),
            TITR_MASS_FO=SIMP(
                statut="f",
                typ=(fonction_sdaster, nappe_sdaster, formule),
                fr=tr("Evolution du titre massique eau vap/eau tot a l'entree [xe]"),
            ),
            PRES_PART=SIMP(
                statut="f", typ="R", fr=tr("Pression partielle air en entree (Pa) [pae]")
            ),
            PRES_PART_FO=SIMP(
                statut="f",
                typ=(fonction_sdaster, nappe_sdaster, formule),
                fr=tr("Evolution de la pression partielle air en entree (Pa) [pae]"),
            ),
        ),
        b_condition_6=BLOC(
            condition="""equal_to("FLUIDE_ENTREE", 6)""",
            regles=(UN_PARMI("TEMP_ENTREE", "TEMP_ENTREE_FO")),
            fr=tr("Air seul"),
            TEMP_ENTREE=SIMP(statut="f", typ="R", fr=tr("Temperature a l'entree (degres C) [te]")),
            TEMP_ENTREE_FO=SIMP(
                statut="f",
                typ=(fonction_sdaster, nappe_sdaster, formule),
                fr=tr("Evolution de la temperature a l'entree (degres C) [te]"),
            ),
        ),
    ),
    LIST_INST=SIMP(statut="f", typ=(listr8_sdaster), fr=tr("Liste des instants de calcul imposes")),
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
            condition="""equal_to("ECOULEMENT", 'GELE')""",
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
            condition="""less_than('FROTTEMENT', 0)""",
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
            condition="""less_than('TRANSFERT_CHAL', 0)""",
            fr=tr("Cas diphasique"),
            XMINCH=SIMP(statut="o", typ="R", fr=tr("Titre massique gazeux min [xminch]")),
            XMAXCH=SIMP(statut="o", typ="R", fr=tr("Titre massique gazeux max [xmaxch]")),
        ),
    ),
    #      CRITERE DE CONVERGENCE
    #      **********************
    CONV_CRITERE=FACT(
        statut="o",
        min=1,
        max=1,
        TEMP_REF=SIMP(
            statut="o",
            typ="R",
            val_min=1.0e-5,
            fr=tr("Temperature de reference pour le calcul du critere"),
        ),
        PRES_REF=SIMP(
            statut="o",
            typ="R",
            val_min=1.0e-5,
            fr=tr("Pression de reference pour le calcul du critere"),
        ),
        CRITERE=SIMP(
            statut="f",
            typ="TXM",
            defaut="TEMP_PRESS",
            into=("TEMP_PRESS", "EXPLICITE", "TEMP", "PRESS"),
            fr=tr("La nature du critere pour la convergence"),
        ),
        b_critere_autre=BLOC(
            condition="""equal_to("CRITERE", 'TEMP_PRESS') or equal_to("CRITERE", 'TEMP') or equal_to("CRITERE", 'PRESS')""",
            fr=tr("Critere de convergence temp_press, temp, ou press"),
            SUBD_NIVEAU=SIMP(
                statut="f",
                typ="I",
                val_min=2,
                defaut=3,
                fr=tr("Nombre maximum de niveau de subdivision d'un pas de temps"),
            ),
            SUBD_PAS_MINI=SIMP(
                statut="f",
                typ="R",
                val_min=0.0,
                fr=tr("Pas de temps en dessous duquel on ne subdivise plus"),
            ),
            NUME_ORDRE_MIN=SIMP(
                statut="f",
                typ="I",
                val_min=-1,
                defaut=-1,
                fr=tr("Numero d'ordre a partir duquel le critere est pris en compte"),
            ),
            PREC_CRIT=SIMP(
                statut="f",
                typ="R",
                val_min=1.0e-2,
                defaut=1.0,
                fr=tr("Valeur du critere pour l'erreur de convergence"),
            ),
        ),
    ),
    #      DONNEES RELATIVES A LA CONVERGENCE NUMERIQUE
    #      ********************************************
    CONVERGENCE_ECREVISSE=FACT(
        statut="f",
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
    ENTETE=SIMP(statut="f", typ="TXM", defaut="Titre du calcul Ecrevisse"),
    IMPRESSION=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
    #      DONNEES POUR STAT_NON_LINE ET THER_NON_LINE
    #      *******************************************
    # copie de stat_non_line.capy des options des mots cles qui nous interessent
    # donnees communes
    CHAM_MATER=SIMP(statut="o", typ=cham_mater),
    TEMP_INIT=SIMP(statut="o", typ="R"),
    CARA_ELEM=SIMP(statut="f", typ=cara_elem),
    # donnees specifiques a stat_non_line
    EXCIT_MECA=FACT(
        statut="o",
        max="**",
        CHARGE=SIMP(statut="o", typ=(char_meca, char_cine_meca)),
        FONC_MULT=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        TYPE_CHARGE=SIMP(
            statut="f", typ="TXM", defaut="FIXE_CSTE", into=("FIXE_CSTE", "SUIV", "DIDI")
        ),
    ),
    CONTACT=SIMP(statut="o", typ=char_contact),
    COMPORTEMENT=C_COMPORTEMENT("MECA_NON_LINE"),
    NEWTON=FACT(
        statut="d",
        REAC_INCR=SIMP(statut="f", typ="I", defaut=1),
        PREDICTION=SIMP(
            statut="f", typ="TXM", into=("DEPL_CALCULE", "TANGENTE", "ELASTIQUE", "EXTRAPOL")
        ),
        MATRICE=SIMP(statut="f", typ="TXM", defaut="TANGENTE", into=("TANGENTE", "ELASTIQUE")),
        PAS_MINI_ELAS=SIMP(statut="f", typ="R", defaut=0.0e0),
        REAC_ITER=SIMP(statut="f", typ="I", defaut=0),
        REAC_ITER_ELAS=SIMP(statut="f", typ="I", defaut=0),
        EVOL_NOLI=SIMP(statut="f", typ=evol_noli),
    ),
    CONVERGENCE=FACT(
        statut="d",
        regles=(PRESENT_ABSENT("RESI_REFE_RELA", "RESI_GLOB_MAXI", "RESI_GLOB_RELA"),),
        b_refe_rela=BLOC(
            condition="""exists("RESI_REFE_RELA")""",
            regles=(
                AU_MOINS_UN(
                    "SIGM_REFE",
                    "EPSI_REFE",
                    "FLUX_THER_REFE",
                    "FLUX_HYD1_REFE",
                    "FLUX_HYD2_REFE",
                    "VARI_REFE",
                ),
            ),
            SIGM_REFE=SIMP(statut="f", typ="R"),
            EPSI_REFE=SIMP(statut="f", typ="R"),
            FLUX_THER_REFE=SIMP(statut="f", typ="R"),
            FLUX_HYD1_REFE=SIMP(statut="f", typ="R"),
            FLUX_HYD2_REFE=SIMP(statut="f", typ="R"),
            VARI_REFE=SIMP(statut="f", typ="R"),
        ),
        RESI_REFE_RELA=SIMP(statut="f", typ="R"),
        RESI_GLOB_MAXI=SIMP(statut="f", typ="R"),
        RESI_GLOB_RELA=SIMP(statut="f", typ="R"),
        ITER_GLOB_MAXI=SIMP(statut="f", typ="I", defaut=10),
        ITER_GLOB_ELAS=SIMP(statut="f", typ="I", defaut=25),
        ARRET=SIMP(statut="f", typ="TXM", defaut="OUI", into=("OUI", "NON")),
    ),
    ENERGIE=FACT(
        statut="f", max=1, CALCUL=SIMP(statut="f", typ="TXM", into=("OUI",), defaut="OUI")
    ),
    # donnees specifiques a ther_lineaire
    EXCIT_THER=FACT(
        statut="o",
        max="**",
        CHARGE=SIMP(statut="o", typ=(char_ther, char_cine_ther)),
        FONC_MULT=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    ),
)
