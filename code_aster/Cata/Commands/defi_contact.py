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

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *

DEFI_CONTACT = OPER(
    nom="DEFI_CONTACT",
    op=30,
    sd_prod=char_contact,
    reentrant="n",
    fr=tr(
        "Définit les zones soumises à des conditions de contact unilatéral avec ou sans frottement"
    ),
    # en        = "Allows the definition of contact surfaces as well as unilateral conditions",
    # MODELE
    MODELE=SIMP(statut="o", typ=modele_sdaster),
    INFO=SIMP(statut="f", typ="I", into=(1, 2)),
    # FORMULATION (UNIQUE PAR OCCURRENCE DE DEFI_CONTACT)
    FORMULATION=SIMP(
        statut="f",
        typ="TXM",
        fr=tr("Choix d'une formulation de contact ou de liaisons unilatérales"),
        defaut="DISCRETE",
        into=("DISCRETE", "CONTINUE", "LIAISON_UNIL"),
    ),
    # PARAMETRE GENERAL : FROTTEMENT
    FROTTEMENT=SIMP(
        statut="f",
        typ="TXM",
        defaut="SANS",
        into=("COULOMB", "SANS"),
        fr=tr("Choix d'un modèle de frottement (uniquement pour les formulations de contact)"),
    ),
    ### PARAMETRES GENERAUX (UNIQUEMENT POUR LE CONTACT MAILLE, NE DEPENDENT PAS DE LA ZONE DE CONTACT)
    b_contact_mail=BLOC(
        condition="""((equal_to("FORMULATION", 'CONTINUE')) or (equal_to("FORMULATION", 'DISCRETE')))""",
        #          ARRET DU CALCUL POUR LE MODE SANS RESOLUTION DU CONTACT
        STOP_INTERP=SIMP(
            statut="f",
            typ="TXM",
            defaut="NON",
            into=("OUI", "NON"),
            fr=tr(
                "Arrête le calcul dès qu'une interpénétration est détectée en mode RESOLUTION='NON'"
            ),
        ),
        #          LISSAGE DES NORMALES PAR MOYENNATION AUX NOEUDS
        LISSAGE=SIMP(
            statut="f",
            typ="TXM",
            defaut="NON",
            into=("OUI", "NON"),
            fr=tr("Lissage des normales par moyennation aux noeuds"),
        ),
        #          VERIFICATION DE L'ORIENTATION ET DE LA COHERENCE DES NORMALES
        VERI_NORM=SIMP(
            statut="f",
            typ="TXM",
            defaut="OUI",
            into=("OUI", "NON"),
            fr=tr("Vérification de l'orientation (sortante) des normales aux surfaces"),
        ),
    ),
    ### PARAMETRES GENERAUX (UNIQUEMENT POUR LE CONTACT, NE DEPENDENT PAS DE LA ZONE DE CONTACT)
    # PARAMETRE GENERAL : BOUCLE DE GEOMETRIE - Cas discret
    b_bouc_geom_maild=BLOC(
        condition="""equal_to("FORMULATION", 'DISCRETE')""",
        ALGO_RESO_GEOM=SIMP(statut="f", typ="TXM", into=("POINT_FIXE",), defaut="POINT_FIXE"),
        REAC_GEOM=SIMP(
            statut="f", typ="TXM", into=("AUTOMATIQUE", "CONTROLE", "SANS"), defaut="AUTOMATIQUE"
        ),
        b_automatique=BLOC(
            condition="""equal_to("REAC_GEOM", 'AUTOMATIQUE') """,
            ITER_GEOM_MAXI=SIMP(statut="f", typ="I", defaut=10),
            RESI_GEOM=SIMP(statut="f", typ="R", defaut=0.01),
        ),
        b_controle=BLOC(
            condition="""equal_to("REAC_GEOM", 'CONTROLE') """,
            NB_ITER_GEOM=SIMP(statut="f", typ="I", defaut=2),
        ),
    ),
    # PARAMETRE GENERAL : BOUCLE DE GEOMETRIE - Cas continu
    b_bouc_geom_mailc=BLOC(
        condition="""equal_to("FORMULATION", 'CONTINUE')""",
        ALGO_RESO_GEOM=SIMP(
            statut="f", typ="TXM", into=("POINT_FIXE", "NEWTON"), defaut="POINT_FIXE"
        ),
        CONT_STAT_ELAS=SIMP(statut="f", typ="I", val_min=0, defaut=0),
        b_algo_reso_geomNE=BLOC(
            condition="""equal_to("ALGO_RESO_GEOM", 'NEWTON')""",
            RESI_GEOM=SIMP(statut="f", typ="R", defaut=0.000001),
        ),
        b_algo_reso_geomPF=BLOC(
            condition="""equal_to("ALGO_RESO_GEOM", 'POINT_FIXE')""",
            REAC_GEOM=SIMP(
                statut="f",
                typ="TXM",
                into=("AUTOMATIQUE", "CONTROLE", "SANS"),
                defaut="AUTOMATIQUE",
            ),
            b_automatique=BLOC(
                condition="""equal_to("REAC_GEOM", 'AUTOMATIQUE') """,
                ITER_GEOM_MAXI=SIMP(statut="f", typ="I", defaut=10),
                RESI_GEOM=SIMP(statut="f", typ="R", defaut=0.01),
            ),
            b_controle=BLOC(
                condition="""equal_to("REAC_GEOM", 'CONTROLE') """,
                NB_ITER_GEOM=SIMP(statut="f", typ="I", defaut=2),
            ),
        ),
    ),
    # PARAMETRE GENERAL : BOUCLE DE CONTACT
    b_bouc_cont_disc=BLOC(
        condition="""equal_to("FORMULATION", 'DISCRETE') """,
        ITER_CONT_MULT=SIMP(statut="f", typ="I", defaut=4),
    ),
    b_bouc_cont_cont=BLOC(
        condition="""equal_to("FORMULATION", 'CONTINUE') """,
        ALGO_RESO_CONT=SIMP(statut="f", typ="TXM", defaut="NEWTON", into=("POINT_FIXE", "NEWTON")),
        RESI_CONT=SIMP(statut="f", typ="R", defaut=-1),
        b_algo_reso_contPF=BLOC(
            condition="""equal_to("ALGO_RESO_CONT", 'POINT_FIXE')""",
            ITER_CONT_TYPE=SIMP(statut="f", typ="TXM", defaut="MAXI", into=("MULT", "MAXI")),
            b_bouc_cont_mult=BLOC(
                condition="""equal_to("ITER_CONT_TYPE", 'MULT')""",
                ITER_CONT_MULT=SIMP(statut="f", typ="I", defaut=4),
            ),
            b_bouc_cont_maxi=BLOC(
                condition="""equal_to("ITER_CONT_TYPE", 'MAXI')""",
                ITER_CONT_MAXI=SIMP(statut="f", typ="I", defaut=30),
            ),
        ),
    ),
    # PARAMETRE GENERAL : BOUCLE DE FROTTEMENT - Cas continu
    b_bouc_frot_cont=BLOC(
        condition="""equal_to("FROTTEMENT", 'COULOMB') and equal_to("FORMULATION", 'CONTINUE') """,
        ALGO_RESO_FROT=SIMP(statut="f", typ="TXM", into=("POINT_FIXE", "NEWTON"), defaut="NEWTON"),
        b_algo_reso_frotPF=BLOC(
            condition="""equal_to("ALGO_RESO_FROT", 'POINT_FIXE')""",
            ITER_FROT_MAXI=SIMP(statut="f", typ="I", defaut=10),
            RESI_FROT=SIMP(statut="f", typ="R", defaut=0.0001),
        ),
        b_algo_reso_frotNE=BLOC(
            condition="""equal_to("ALGO_RESO_FROT", 'NEWTON')""",
            RESI_FROT=SIMP(statut="f", typ="R", defaut=0.0001),
        ),
    ),
    # PARAMETRES GENERAUX : METHODES DISCRETES
    b_para_discret=BLOC(
        condition="""equal_to("FORMULATION", 'DISCRETE') """,
        #                                 ## METHODES DE DUALISATION ##
        STOP_SINGULIER=SIMP(
            statut="f",
            typ="TXM",
            fr=tr("Tient compte de la singularité de la matrice de contact"),
            defaut="OUI",
            into=("OUI", "NON"),
        ),
        NB_RESOL=SIMP(
            statut="f",
            typ="I",
            fr=tr("Nombre de résolutions simultanées pour la construction du complément de Schur"),
            defaut=10,
        ),
        #                                 ## GCP ##
        RESI_ABSO=SIMP(
            statut="f",
            typ="R",
            fr=tr("Critère de convergence (niveau d'interpénétration autorisé pour 'GCP')"),
        ),
        ITER_GCP_MAXI=SIMP(
            statut="f", typ="I", fr=tr("Nombre d'itérations maximal ('GCP')"), defaut=0
        ),
        RECH_LINEAIRE=SIMP(
            statut="f",
            typ="TXM",
            fr=tr("Autorisation de sortie du domaine admissible lors de la recherche linéaire"),
            defaut="ADMISSIBLE",
            into=("ADMISSIBLE", "NON_ADMISSIBLE"),
        ),
        PRE_COND=SIMP(
            statut="f",
            typ="TXM",
            fr=tr("Choix d'un préconditionneur (accélère la convergence de 'GCP')"),
            defaut="SANS",
            into=("DIRICHLET", "SANS"),
        ),
        b_dirichlet=BLOC(
            condition="""equal_to("PRE_COND", 'DIRICHLET')""",
            COEF_RESI=SIMP(
                statut="f",
                typ="R",
                fr=tr("Activation du préconditionneur quand le résidu a été divisé par COEF_RESI"),
                defaut=-1.0,
            ),
            ITER_PRE_MAXI=SIMP(
                statut="f",
                typ="I",
                fr=tr("Nombre d'itérations maximal pour le préconditionneur ('GCP')"),
                defaut=0,
            ),
        ),
    ),  # fin b_para_discret
    ## AFFECTATIONS (ZONES PAR ZONES)
    # AFFECTATION - CAS LIAISON_UNILATERALE
    b_affe_unil=BLOC(
        condition="""equal_to("FORMULATION", 'LIAISON_UNIL')""",
        ZONE=FACT(
            statut="o",
            max="**",
            # -- Liaison unilatérale
            regles=(UN_PARMI("GROUP_MA", "GROUP_NO"),),
            GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
            GROUP_NO=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
            #
            NOM_CMP=SIMP(statut="o", typ="TXM", max="**"),
            COEF_IMPO=SIMP(statut="o", typ=(fonction_sdaster, nappe_sdaster, formule)),
            COEF_MULT=SIMP(statut="o", typ=(fonction_sdaster, nappe_sdaster, formule), max="**"),
            # -- Incompatibilité avec CL
            SANS_GROUP_NO=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
            # --- Résolution
            ALGO_CONT=SIMP(
                statut="f", typ="TXM", defaut="CONTRAINTE", into=("CONTRAINTE", "PENALISATION")
            ),
            b_cont_pen=BLOC(
                condition="""equal_to("ALGO_CONT", 'PENALISATION') """,
                fr=tr("Paramètre de la méthode pénalisée"),
                COEF_PENA=SIMP(statut="o", typ="R"),
            ),
        ),  # fin mot-clé facteur ZONE
    ),  # fin b_affe_unil
    # AFFECTATION - CAS DISCRET et FROTTEMENT='COULOMB'
    b_affe_discret=BLOC(
        condition="""equal_to("FORMULATION", 'DISCRETE') and equal_to("FROTTEMENT", 'COULOMB')""",
        ZONE=FACT(
            statut="o",
            max="**",
            # --- Appariement
            APPARIEMENT=SIMP(
                statut="f", typ="TXM", defaut="MAIT_ESCL", into=("NODAL", "MAIT_ESCL")
            ),
            #
            GROUP_MA_MAIT=SIMP(statut="f", typ=grma, max=1),
            GROUP_MA_ESCL=SIMP(statut="f", typ=grma, max=1),
            #
            NORMALE=SIMP(statut="f", typ="TXM", defaut="MAIT", into=("MAIT", "MAIT_ESCL", "ESCL")),
            #
            VECT_MAIT=SIMP(statut="f", typ="TXM", defaut="AUTO", into=("AUTO", "FIXE", "VECT_Y")),
            b_nmait_fixe=BLOC(
                condition="""equal_to("VECT_MAIT", 'FIXE')""",
                MAIT_FIXE=SIMP(statut="o", typ="R", min=3, max=3),
            ),
            b_nmait_vecty=BLOC(
                condition="""equal_to("VECT_MAIT", 'VECT_Y')""",
                MAIT_VECT_Y=SIMP(statut="o", typ="R", min=3, max=3),
            ),
            #
            VECT_ESCL=SIMP(statut="f", typ="TXM", defaut="AUTO", into=("AUTO", "FIXE", "VECT_Y")),
            b_nescl_fixe=BLOC(
                condition="""equal_to("VECT_ESCL", 'FIXE')""",
                ESCL_FIXE=SIMP(statut="o", typ="R", min=3, max=3),
            ),
            b_nescl_vecty=BLOC(
                condition="""equal_to("VECT_ESCL", 'VECT_Y')""",
                ESCL_VECT_Y=SIMP(statut="o", typ="R", min=3, max=3),
            ),
            #
            TYPE_PROJECTION=SIMP(
                statut="f", typ="TXM", defaut="ORTHOGONALE", into=("ORTHOGONALE", "FIXE")
            ),
            b_appa_fixe=BLOC(
                condition="""equal_to("TYPE_PROJECTION", 'FIXE')""",
                DIRE_APPA=SIMP(statut="f", typ="R", min=3, max=3),
            ),
            #
            DIST_POUTRE=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
            DIST_COQUE=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
            b_cara=BLOC(
                condition="""equal_to("DIST_POUTRE", 'OUI') or equal_to("DIST_COQUE", 'OUI')""",
                CARA_ELEM=SIMP(statut="o", typ=(cara_elem)),
            ),
            DIST_MAIT=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
            DIST_ESCL=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
            #
            DIST_APPA=SIMP(statut="f", typ="R", defaut=-1.0),
            TOLE_PROJ_EXT=SIMP(statut="f", typ="R", defaut=0.50),
            # --- Incompatibilité avec CL
            SANS_GROUP_NO=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
            SANS_GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
            # --- Mode sans calcul
            RESOLUTION=SIMP(statut="f", typ="TXM", defaut="OUI", into=("OUI", "NON")),
            b_verif=BLOC(
                condition="""equal_to("RESOLUTION", 'NON') """,
                TOLE_INTERP=SIMP(statut="f", typ="R", defaut=0.0, val_min=0.0),
            ),
            # --- Résolution
            ALGO_CONT=SIMP(
                statut="f",
                typ="TXM",
                defaut="CONTRAINTE",
                into=("CONTRAINTE", "PENALISATION", "GCP"),
            ),
            b_active=BLOC(
                condition="""equal_to("ALGO_CONT", 'CONTRAINTE') """,
                fr=tr("Paramètres de la méthode des contraintes actives (contact uniquement)"),
                GLISSIERE=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
                b_glissiere=BLOC(
                    condition="""equal_to("GLISSIERE", 'OUI') """,
                    ALARME_JEU=SIMP(statut="f", typ="R", defaut=0.0),
                ),
            ),
            #
            b_penal_contact=BLOC(
                condition="""equal_to("ALGO_CONT", 'PENALISATION') """,
                fr=tr("Paramètres de la méthode pénalisée (contact)"),
                E_N=SIMP(statut="o", typ="R"),
            ),
            #
            COULOMB=SIMP(statut="o", typ="R"),
            COEF_MATR_FROT=SIMP(statut="f", typ="R", defaut=0.0e0),
            ALGO_FROT=SIMP(statut="f", typ="TXM", defaut="PENALISATION", into=("PENALISATION",)),
            E_T=SIMP(statut="o", typ="R"),
        ),  # fin mot-clé facteur ZONE
    ),  # fin bloc b_affe_discret
    b_affe_discret_frot=BLOC(
        condition="""equal_to("FORMULATION", 'DISCRETE') and not equal_to("FROTTEMENT", 'COULOMB')""",
        ZONE=FACT(
            statut="o",
            max="**",
            # --- Appariement
            APPARIEMENT=SIMP(
                statut="f", typ="TXM", defaut="MAIT_ESCL", into=("NODAL", "MAIT_ESCL")
            ),
            #
            GROUP_MA_MAIT=SIMP(statut="f", typ=grma, max=1),
            GROUP_MA_ESCL=SIMP(statut="f", typ=grma, max=1),
            #
            NORMALE=SIMP(statut="f", typ="TXM", defaut="MAIT", into=("MAIT", "MAIT_ESCL", "ESCL")),
            #
            VECT_MAIT=SIMP(statut="f", typ="TXM", defaut="AUTO", into=("AUTO", "FIXE", "VECT_Y")),
            b_nmait_fixe=BLOC(
                condition="""equal_to("VECT_MAIT", 'FIXE')""",
                MAIT_FIXE=SIMP(statut="o", typ="R", min=3, max=3),
            ),
            b_nmait_vecty=BLOC(
                condition="""equal_to("VECT_MAIT", 'VECT_Y')""",
                MAIT_VECT_Y=SIMP(statut="o", typ="R", min=3, max=3),
            ),
            #
            VECT_ESCL=SIMP(statut="f", typ="TXM", defaut="AUTO", into=("AUTO", "FIXE", "VECT_Y")),
            b_nescl_fixe=BLOC(
                condition="""equal_to("VECT_ESCL", 'FIXE')""",
                ESCL_FIXE=SIMP(statut="o", typ="R", min=3, max=3),
            ),
            b_nescl_vecty=BLOC(
                condition="""equal_to("VECT_ESCL", 'VECT_Y')""",
                ESCL_VECT_Y=SIMP(statut="o", typ="R", min=3, max=3),
            ),
            #
            TYPE_PROJECTION=SIMP(
                statut="f", typ="TXM", defaut="ORTHOGONALE", into=("ORTHOGONALE", "FIXE")
            ),
            b_appa_fixe=BLOC(
                condition="""equal_to("TYPE_PROJECTION", 'FIXE')""",
                DIRE_APPA=SIMP(statut="f", typ="R", min=3, max=3),
            ),
            #
            DIST_POUTRE=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
            DIST_COQUE=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
            b_cara=BLOC(
                condition="""equal_to("DIST_POUTRE", 'OUI') or equal_to("DIST_COQUE", 'OUI')""",
                CARA_ELEM=SIMP(statut="o", typ=(cara_elem)),
            ),
            DIST_MAIT=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
            DIST_ESCL=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
            #
            DIST_APPA=SIMP(statut="f", typ="R", defaut=-1.0),
            TOLE_PROJ_EXT=SIMP(statut="f", typ="R", defaut=0.50),
            # --- Incompatibilité avec CL
            SANS_GROUP_NO=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
            SANS_GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
            # --- Mode sans calcul
            RESOLUTION=SIMP(statut="f", typ="TXM", defaut="OUI", into=("OUI", "NON")),
            b_verif=BLOC(
                condition="""equal_to("RESOLUTION", 'NON') """,
                TOLE_INTERP=SIMP(statut="f", typ="R", defaut=0.0, val_min=0.0),
            ),
            # --- Résolution
            ALGO_CONT=SIMP(
                statut="f",
                typ="TXM",
                defaut="CONTRAINTE",
                into=("CONTRAINTE", "PENALISATION", "GCP"),
            ),
            b_active=BLOC(
                condition="""equal_to("ALGO_CONT", 'CONTRAINTE') """,
                fr=tr("Paramètres de la méthode des contraintes actives (contact uniquement)"),
                GLISSIERE=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
                b_glissiere=BLOC(
                    condition="""equal_to("GLISSIERE", 'OUI') """,
                    ALARME_JEU=SIMP(statut="f", typ="R", defaut=0.0),
                ),
            ),
            #
            b_penal_contact=BLOC(
                condition="""equal_to("ALGO_CONT", 'PENALISATION') """,
                fr=tr("Paramètres de la méthode pénalisée (contact)"),
                E_N=SIMP(statut="o", typ="R"),
            ),
        ),  # fin mot-clé facteur ZONE
    ),  # fin bloc b_affe_discret_frot
    # AFFECTATION - CAS CONTINUE et FROTTEMENT = COULOMB
    b_affe_continue_frot=BLOC(
        condition="""equal_to("FORMULATION", 'CONTINUE') and equal_to("FROTTEMENT", 'COULOMB')""",
        ZONE=FACT(
            statut="o",
            max="**",
            # --- Method for contact
            ALGO_CONT=SIMP(
                statut="f", typ="TXM", defaut="STANDARD", into=("STANDARD", "PENALISATION", "LAC")
            ),
            ALGO_FROT=SIMP(
                statut="f", typ="TXM", defaut="STANDARD", into=("STANDARD", "PENALISATION")
            ),
            b_frot_std=BLOC(
                condition="""equal_to("ALGO_FROT", 'STANDARD') """,
                fr=tr("Paramètres de la formulation Lagrangienne"),
            ),
            b_contpen_frotstd=BLOC(
                condition="""(equal_to("ALGO_CONT", 'PENALISATION') and equal_to("ALGO_FROT", 'STANDARD') )""",
                fr=tr("Adaptation pour STANDARD et PENALISATION"),
                ADAPTATION=SIMP(
                    statut="f",
                    typ="TXM",
                    defaut="CYCLAGE",
                    into=("ADAPT_COEF", "CYCLAGE", "TOUT", "NON"),
                ),
                COEF_FROT=SIMP(statut="f", typ="R", defaut=100.0e0),
                b_cont_pena_noadapt=BLOC(
                    condition="""(equal_to("ADAPTATION", 'NON') or equal_to("ADAPTATION", 'CYCLAGE') or not exists("ADAPTATION")) """,
                    fr=tr("Paramètres de la méthode pénalisée"),
                    COEF_PENA_CONT=SIMP(statut="o", typ="R"),
                ),
                b_cont_pena_adapt=BLOC(
                    condition="""(equal_to("ADAPTATION", 'ADAPT_COEF') or equal_to("ADAPTATION", 'TOUT'))  """,
                    fr=tr("Paramètres de la méthode pénalisée"),
                    PENE_MAXI=SIMP(statut="f", typ="R"),
                    COEF_PENA_CONT=SIMP(statut="f", typ="R", defaut=100.0e0),
                ),
            ),
            b_contstd_frotstd=BLOC(
                condition="""(equal_to("ALGO_CONT", 'STANDARD') and equal_to("ALGO_FROT",'STANDARD' ) )""",
                fr=tr("Paramètres de la formulation Lagrangienne"),
                ADAPTATION=SIMP(
                    statut="f",
                    typ="TXM",
                    defaut="CYCLAGE",
                    into=("ADAPT_COEF", "CYCLAGE", "TOUT", "NON"),
                ),
                COEF_CONT=SIMP(statut="f", typ="R", defaut=100.0e0),
                COEF_FROT=SIMP(statut="f", typ="R", defaut=100.0e0),
            ),
            b_contpen_frotpen=BLOC(
                condition="""(equal_to("ALGO_CONT", 'PENALISATION') and equal_to("ALGO_FROT", 'PENALISATION'))""",
                fr=tr("Adaptation pour STANDARD et PENALISATION"),
                ADAPTATION=SIMP(
                    statut="f",
                    typ="TXM",
                    defaut="CYCLAGE",
                    into=("ADAPT_COEF", "CYCLAGE", "TOUT", "NON"),
                ),
                b_cont_pena_noadapt2=BLOC(
                    condition="""(equal_to("ADAPTATION", 'NON') or equal_to("ADAPTATION", 'CYCLAGE') or not exists("ADAPTATION")) """,
                    fr=tr("Paramètres de la méthode pénalisée"),
                    COEF_PENA_CONT=SIMP(statut="o", typ="R"),
                    COEF_PENA_FROT=SIMP(statut="o", typ="R"),
                ),
                b_cont_pena_adapt2=BLOC(
                    condition="""(equal_to("ADAPTATION", 'ADAPT_COEF') or equal_to("ADAPTATION", 'TOUT'))  """,
                    fr=tr("Paramètres de la méthode pénalisée"),
                    PENE_MAXI=SIMP(statut="f", typ="R"),
                    COEF_PENA_CONT=SIMP(statut="f", typ="R", defaut=100.0e0),
                    GLIS_MAXI=SIMP(statut="f", typ="R"),
                    COEF_PENA_FROT=SIMP(statut="f", typ="R", defaut=100.0e5),
                ),
            ),
            b_contstd_frotpen=BLOC(
                condition="""(equal_to("ALGO_CONT", 'STANDARD') and equal_to("ALGO_FROT", 'PENALISATION'))""",
                fr=tr("Paramètres de la formulation Lagrangienne"),
                ADAPTATION=SIMP(
                    statut="f",
                    typ="TXM",
                    defaut="CYCLAGE",
                    into=("ADAPT_COEF", "CYCLAGE", "TOUT", "NON"),
                ),
                COEF_CONT=SIMP(statut="f", typ="R", defaut=100.0e0),
                b_frot_pena_noadapt=BLOC(
                    condition="""(equal_to("ADAPTATION", 'NON') or equal_to("ADAPTATION", 'CYCLAGE'))""",
                    fr=tr("Paramètres de la méthode pénalisée"),
                    COEF_PENA_FROT=SIMP(statut="o", typ="R"),
                ),
                b_frot_pena_adapt=BLOC(
                    condition="""(equal_to("ADAPTATION", 'ADAPT_COEF') or equal_to("ADAPTATION", 'TOUT'))""",
                    fr=tr("Paramètres de la méthode pénalisée"),
                    GLIS_MAXI=SIMP(statut="f", typ="R"),
                    COEF_PENA_FROT=SIMP(statut="f", typ="R", defaut=100.0e5),
                ),
            ),
            # --- Pairing options (for standard slave/master)
            b_pair_std=BLOC(
                condition="""equal_to("ALGO_CONT", 'STANDARD') or equal_to("ALGO_CONT", 'PENALISATION')""",
                fr=tr("Options d'appariement standard "),
                APPARIEMENT=SIMP(statut="f", typ="TXM", defaut="MAIT_ESCL", into=("MAIT_ESCL",)),
                GROUP_MA_MAIT=SIMP(statut="f", typ=grma, max=1),
                GROUP_MA_ESCL=SIMP(statut="f", typ=grma, max=1),
                # -------------------- Choice of normals --------------------------------------------------------------------------------------------------
                NORMALE=SIMP(
                    statut="f", typ="TXM", defaut="MAIT", into=("MAIT", "MAIT_ESCL", "ESCL")
                ),
                VECT_MAIT=SIMP(
                    statut="f", typ="TXM", defaut="AUTO", into=("AUTO", "FIXE", "VECT_Y")
                ),
                b_nmait_fixe=BLOC(
                    condition="""equal_to("VECT_MAIT", 'FIXE')""",
                    MAIT_FIXE=SIMP(statut="o", typ="R", min=3, max=3),
                ),
                b_nmait_vecty=BLOC(
                    condition="""equal_to("VECT_MAIT", 'VECT_Y')""",
                    MAIT_VECT_Y=SIMP(statut="o", typ="R", min=3, max=3),
                ),
                VECT_ESCL=SIMP(
                    statut="f", typ="TXM", defaut="AUTO", into=("AUTO", "FIXE", "VECT_Y")
                ),
                b_nescl_fixe=BLOC(
                    condition="""equal_to("VECT_ESCL", 'FIXE')""",
                    ESCL_FIXE=SIMP(statut="o", typ="R", min=3, max=3),
                ),
                b_nescl_vecty=BLOC(
                    condition="""equal_to("VECT_ESCL", 'VECT_Y')""",
                    ESCL_VECT_Y=SIMP(statut="o", typ="R", min=3, max=3),
                ),
                # -------------------- Options for pairing  -----------------------------------------------------------------------------------------------
                TYPE_PROJECTION=SIMP(
                    statut="f", typ="TXM", defaut="ORTHOGONALE", into=("ORTHOGONALE", "FIXE")
                ),
                b_appa_fixe=BLOC(
                    condition="""equal_to("TYPE_PROJECTION", 'FIXE')""",
                    DIRE_APPA=SIMP(statut="f", typ="R", min=3, max=3),
                ),
                DIST_APPA=SIMP(statut="f", typ="R", defaut=-1.0),
                TOLE_PROJ_EXT=SIMP(statut="f", typ="R", defaut=0.50),
                # -------------------- Add suppl. gaps      -----------------------------------------------------------------------------------------------
                DIST_POUTRE=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
                DIST_COQUE=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
                b_cara=BLOC(
                    condition="""equal_to("DIST_POUTRE", 'OUI') or equal_to("DIST_COQUE", 'OUI')""",
                    CARA_ELEM=SIMP(statut="o", typ=(cara_elem)),
                ),
                DIST_MAIT=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
                DIST_ESCL=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
                # -------------------- Managing boundary conditions with contact --------------------------------------------------------------------------
                SANS_GROUP_NO=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
                SANS_GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
                # -------------------- Specific collocation integration rules -----------------------------------------------------------------------------
                INTEGRATION=SIMP(
                    statut="f",
                    typ="TXM",
                    defaut="AUTO",
                    into=("AUTO", "GAUSS", "SIMPSON", "NCOTES"),
                ),
                b_gauss=BLOC(
                    condition="""equal_to("INTEGRATION", 'GAUSS') """,
                    fr=tr("Degré du polynôme de Legendre donnant les points de Gauss"),
                    ORDRE_INT=SIMP(statut="f", typ="I", defaut=3, val_min=1, val_max=6),
                ),
                b_simpson=BLOC(
                    condition="""equal_to("INTEGRATION", 'SIMPSON') """,
                    fr=tr("Nombre de subdivisions du domaine"),
                    ORDRE_INT=SIMP(statut="f", typ="I", defaut=1, val_min=1, val_max=4),
                ),
                b_ncotes=BLOC(
                    condition="""equal_to("INTEGRATION", 'NCOTES') """,
                    fr=tr("Degré du polynôme interpolateur"),
                    ORDRE_INT=SIMP(statut="f", typ="I", defaut=3, val_min=3, val_max=8),
                ),
                # -------------------- No contact solving (only pairing)  -----------------------------------------------------------------------------------
                RESOLUTION=SIMP(statut="f", typ="TXM", defaut="OUI", into=("OUI", "NON")),
                b_verif=BLOC(
                    condition="""equal_to("RESOLUTION", 'NON') """,
                    TOLE_INTERP=SIMP(statut="f", typ="R", defaut=0.0),
                ),
                # -------------------- Special options  -----------------------------------------------------------------------------------------------------
                CONTACT_INIT=SIMP(
                    statut="f",
                    typ="TXM",
                    defaut="INTERPENETRE",
                    into=("OUI", "INTERPENETRE", "NON"),
                ),
                b_cont_init_inte=BLOC(
                    condition="""equal_to("CONTACT_INIT", 'INTERPENETRE') """,
                    DIST_MAXI=SIMP(statut="f", typ="R", max=1),
                ),
                GLISSIERE=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
            ),  # fin b_pair_std
            # --- Pairing options (for segment to segment contact)
            b_pair_mortar=BLOC(
                condition="""equal_to("ALGO_CONT", 'LAC') """,
                fr=tr("Options d'appariement standard "),
                APPARIEMENT=SIMP(statut="f", typ="TXM", defaut="MORTAR", into=("MORTAR",)),
                TYPE_APPA=SIMP(statut="f", typ="TXM", defaut="RAPIDE", into=("ROBUSTE", "RAPIDE")),
                RESI_APPA=SIMP(statut="f", typ="R", defaut=1.0e-8),
                TYPE_JACOBIEN=SIMP(
                    statut="f", typ="TXM", defaut="INITIAL", into=("INITIAL", "ACTUALISE")
                ),
                INTEGRATION=SIMP(statut="f", typ="TXM", defaut="AUTO", into=("AUTO",)),
                GROUP_MA_MAIT=SIMP(statut="o", typ=grma, max=1),
                GROUP_MA_ESCL=SIMP(statut="o", typ=grma, max=1),
                CONTACT_INIT=SIMP(
                    statut="f",
                    typ="TXM",
                    defaut="INTERPENETRE",
                    into=("OUI", "INTERPENETRE", "NON"),
                ),
                b_cont_init_inte=BLOC(
                    condition="""equal_to("CONTACT_INIT", 'INTERPENETRE') """,
                    DIST_MAXI=SIMP(statut="f", typ="R", max=1),
                ),
            ),  # fin b_pair_morta
            # --- Method for friction
            COULOMB=SIMP(statut="o", typ="R"),
            GRAND_GLIS=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
            SEUIL_INIT=SIMP(statut="f", typ="R"),
            SANS_GROUP_NO_FR=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
            b_sans_group_no_frot=BLOC(
                condition="""exists("SANS_GROUP_NO_FR")""",
                fr=tr("Direction de frottement à exclure (uniquement dans le cas 3D)"),
                DIRE_EXCL_FROT=SIMP(statut="f", typ="R", min=3, max=3),
            ),
        ),  # fin mot-clé facteur ZONE
    ),  # fin bloc b_affe_continue_frot
    # AFFECTATION - CAS CONTINUE et FROTTEMENT != COULOMB
    b_affe_continue=BLOC(
        condition="""equal_to("FORMULATION", 'CONTINUE') and not equal_to("FROTTEMENT", 'COULOMB')""",
        ZONE=FACT(
            statut="o",
            max="**",
            # --- Method for contact
            ALGO_CONT=SIMP(
                statut="f", typ="TXM", defaut="STANDARD", into=("STANDARD", "PENALISATION", "LAC")
            ),
            b_cont_std=BLOC(
                condition="""equal_to("ALGO_CONT", 'STANDARD') """,
                fr=tr("Paramètres de la formulation Lagrangienne"),
                COEF_CONT=SIMP(statut="f", typ="R", defaut=100.0e0),
                ADAPTATION=SIMP(
                    statut="f",
                    typ="TXM",
                    defaut="CYCLAGE",
                    into=("ADAPT_COEF", "CYCLAGE", "TOUT", "NON"),
                ),
            ),
            b_cont_pen=BLOC(
                condition="""equal_to("ALGO_CONT", 'PENALISATION')  """,
                fr=tr("Paramètres de la méthode pénalisée"),
                # COEF_PENA_CONT  =SIMP(statut='o',typ='R'),
                # COEF_PENA_CONT  =SIMP(statut='f',typ='R',defaut=100.0E+0),
                ADAPTATION=SIMP(
                    statut="f",
                    typ="TXM",
                    defaut="CYCLAGE",
                    into=("ADAPT_COEF", "CYCLAGE", "TOUT", "NON"),
                ),
                # PENE_MAXI  =SIMP(statut='f',typ='R'),
                b_adapt_penal=BLOC(
                    condition="""(equal_to("ADAPTATION", 'ADAPT_COEF')) or (equal_to("ADAPTATION", 'TOUT'))""",
                    PENE_MAXI=SIMP(statut="f", typ="R"),
                    COEF_PENA_CONT=SIMP(statut="f", typ="R", defaut=100.0e0),
                ),
                b_noadapt_penal=BLOC(
                    condition="""(equal_to("ADAPTATION", 'CYCLAGE')) or (equal_to("ADAPTATION", 'NON')) or (not exists("ADAPTATION"))""",
                    # PENE_MAXI  =SIMP(statut='f',typ='R'),
                    COEF_PENA_CONT=SIMP(statut="o", typ="R"),
                ),
            ),
            # --- Pairing options (for standard slave/master)
            b_pair_std=BLOC(
                condition="""equal_to("ALGO_CONT", 'STANDARD') or equal_to("ALGO_CONT", 'PENALISATION') """,
                fr=tr("Options d'appariement standard "),
                APPARIEMENT=SIMP(statut="f", typ="TXM", defaut="MAIT_ESCL", into=("MAIT_ESCL",)),
                GROUP_MA_MAIT=SIMP(statut="f", typ=grma, max=1),
                GROUP_MA_ESCL=SIMP(statut="f", typ=grma, max=1),
                # -------------------- Choice of normals --------------------------------------------------------------------------------------------------
                NORMALE=SIMP(
                    statut="f", typ="TXM", defaut="MAIT", into=("MAIT", "MAIT_ESCL", "ESCL")
                ),
                VECT_MAIT=SIMP(
                    statut="f", typ="TXM", defaut="AUTO", into=("AUTO", "FIXE", "VECT_Y")
                ),
                b_nmait_fixe=BLOC(
                    condition="""equal_to("VECT_MAIT", 'FIXE')""",
                    MAIT_FIXE=SIMP(statut="o", typ="R", min=3, max=3),
                ),
                b_nmait_vecty=BLOC(
                    condition="""equal_to("VECT_MAIT", 'VECT_Y')""",
                    MAIT_VECT_Y=SIMP(statut="o", typ="R", min=3, max=3),
                ),
                VECT_ESCL=SIMP(
                    statut="f", typ="TXM", defaut="AUTO", into=("AUTO", "FIXE", "VECT_Y")
                ),
                b_nescl_fixe=BLOC(
                    condition="""equal_to("VECT_ESCL", 'FIXE')""",
                    ESCL_FIXE=SIMP(statut="o", typ="R", min=3, max=3),
                ),
                b_nescl_vecty=BLOC(
                    condition="""equal_to("VECT_ESCL", 'VECT_Y')""",
                    ESCL_VECT_Y=SIMP(statut="o", typ="R", min=3, max=3),
                ),
                # -------------------- Options for pairing  -----------------------------------------------------------------------------------------------
                TYPE_PROJECTION=SIMP(
                    statut="f", typ="TXM", defaut="ORTHOGONALE", into=("ORTHOGONALE", "FIXE")
                ),
                b_appa_fixe=BLOC(
                    condition="""equal_to("TYPE_PROJECTION", 'FIXE')""",
                    DIRE_APPA=SIMP(statut="f", typ="R", min=3, max=3),
                ),
                DIST_APPA=SIMP(statut="f", typ="R", defaut=-1.0),
                TOLE_PROJ_EXT=SIMP(statut="f", typ="R", defaut=0.50),
                # -------------------- Add suppl. gaps      -----------------------------------------------------------------------------------------------
                DIST_POUTRE=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
                DIST_COQUE=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
                b_cara=BLOC(
                    condition="""equal_to("DIST_POUTRE", 'OUI') or equal_to("DIST_COQUE", 'OUI')""",
                    CARA_ELEM=SIMP(statut="o", typ=(cara_elem)),
                ),
                DIST_MAIT=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
                DIST_ESCL=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
                # -------------------- Managing boundary conditions with contact --------------------------------------------------------------------------
                SANS_GROUP_NO=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
                SANS_GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
                # -------------------- Specific collocation integration rules -----------------------------------------------------------------------------
                INTEGRATION=SIMP(
                    statut="f",
                    typ="TXM",
                    defaut="AUTO",
                    into=("AUTO", "GAUSS", "SIMPSON", "NCOTES"),
                ),
                b_gauss=BLOC(
                    condition="""equal_to("INTEGRATION", 'GAUSS') """,
                    fr=tr("Degré du polynôme de Legendre donnant les points de Gauss"),
                    ORDRE_INT=SIMP(statut="f", typ="I", defaut=3, val_min=1, val_max=6),
                ),
                b_simpson=BLOC(
                    condition="""equal_to("INTEGRATION", 'SIMPSON') """,
                    fr=tr("Nombre de subdivisions du domaine"),
                    ORDRE_INT=SIMP(statut="f", typ="I", defaut=1, val_min=1, val_max=4),
                ),
                b_ncotes=BLOC(
                    condition="""equal_to("INTEGRATION", 'NCOTES') """,
                    fr=tr("Degré du polynôme interpolateur"),
                    ORDRE_INT=SIMP(statut="f", typ="I", defaut=3, val_min=3, val_max=8),
                ),
                # -------------------- No contact solving (only pairing)  -----------------------------------------------------------------------------------
                RESOLUTION=SIMP(statut="f", typ="TXM", defaut="OUI", into=("OUI", "NON")),
                b_verif=BLOC(
                    condition="""equal_to("RESOLUTION", 'NON') """,
                    TOLE_INTERP=SIMP(statut="f", typ="R", defaut=0.0),
                ),
                # -------------------- Special options  -----------------------------------------------------------------------------------------------------
                CONTACT_INIT=SIMP(
                    statut="f",
                    typ="TXM",
                    defaut="INTERPENETRE",
                    into=("OUI", "INTERPENETRE", "NON"),
                ),
                b_cont_init_inte=BLOC(
                    condition="""equal_to("CONTACT_INIT", 'INTERPENETRE') """,
                    DIST_MAXI=SIMP(statut="f", typ="R", max=1),
                ),
                GLISSIERE=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
            ),  # fin b_pair_std
            # --- Pairing options (for segment to segment contact)
            b_pair_mortar=BLOC(
                condition="""equal_to("ALGO_CONT", 'LAC') """,
                fr=tr("Options d'appariement standard "),
                APPARIEMENT=SIMP(statut="f", typ="TXM", defaut="MORTAR", into=("MORTAR",)),
                TYPE_APPA=SIMP(statut="f", typ="TXM", defaut="RAPIDE", into=("ROBUSTE", "RAPIDE")),
                RESI_APPA=SIMP(statut="f", typ="R", defaut=1.0e-8),
                TYPE_JACOBIEN=SIMP(
                    statut="f", typ="TXM", defaut="INITIAL", into=("INITIAL", "ACTUALISE")
                ),
                INTEGRATION=SIMP(statut="f", typ="TXM", defaut="AUTO", into=("AUTO",)),
                GROUP_MA_MAIT=SIMP(statut="o", typ=grma, max=1),
                GROUP_MA_ESCL=SIMP(statut="o", typ=grma, max=1),
                CONTACT_INIT=SIMP(
                    statut="f",
                    typ="TXM",
                    defaut="INTERPENETRE",
                    into=("OUI", "INTERPENETRE", "NON"),
                ),
            ),
        ),  # fin mot-clé facteur ZONE
    ),  # fin bloc b_affe_continue
)  # fin OPER
