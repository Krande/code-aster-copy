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

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *


# quelques remarques:
# - ce catalogue est un exemple de ce à quoi il pourrait ressembler à la fin du chantier en 2024
# - on ne passera plus par un op fortran. Pur c++/python
# - il faudra créer une nouvelle classe (et sd) pour ce contact
# - il faudra ajouter le maximum de fonctions d'intérogation à cette classe (pour rajouter des tests facilement)
# - il faut regarder ce que cela donne dans asterStudy
# - est-ce que l'on veut rajouter des choses qui ne sont possibles actuellement ?

DEFI_CONT=OPER(nom = "DEFI_CONT", op=None, sd_prod   = char_contact, reentrant = 'n',
                fr = tr("Définit les zones soumises à des conditions de contact avec ou sans frottement"),
                #en        = "Allows the definition of contact surfaces",

# ----- PARAMETRES GENERAUX ( NE DEPENDENT PAS DE LA ZONE DE CONTACT)

        MODELE          =SIMP(statut='o',typ=modele_sdaster,),
        INFO            =SIMP(statut='f',typ='I',into=(1,2),defaut=1),

# PARAMETRE GENERAL : FROTTEMENT (LE CHOIX DU TYPE EST FAIT PAR ZONE APRES)
        FROTTEMENT      =SIMP(statut='f', typ='TXM', defaut="NON", into=("NON", "OUI"),
                               fr=tr("Activation du frottement" ),),

# LISSAGE DES NORMALES PAR MOYENNATION AUX NOEUDS
        LISSAGE         = SIMP(statut='f', typ='TXM', defaut="NON", into=("OUI","NON"),
                                fr=tr("Lissage des normales par moyennation aux noeuds"),),


# ----- définition mot-clé facteur sans frottement
        b_zone_nofric     = BLOC(condition = """equal_to("FROTTEMENT", 'NON') """,
                ZONE = FACT(statut='o', max='**',
# VERIFICATION DE L'ORIENTATION ET DE LA COHERENCE DES NORMALES
                VERI_NORM       =SIMP(statut='f', typ='TXM', defaut="OUI", into=("OUI","NON"),
                                 fr=tr("Vérification de l'orientation (sortante) des normales aux surfaces"),),
# Method for contact
                ALGO_CONT   = SIMP(statut='f',typ='TXM',defaut="LAGRANGIEN",
                                        into=("LAGRANGIEN","NITSCHE","PENALISATION",),),

                b_algo_cont = BLOC(condition = """is_in("ALGO_CONT", ('LAGRANGIEN', 'NITSCHE'))""",
                                   VARIANTE  =SIMP(statut='f', typ='TXM', defaut="ROBUSTE",
                                                into=("RAPIDE","ROBUSTE"),
                                                fr=tr("Choix de la variante des formulations du contact" ),),
                                   b_vari_syme = BLOC(condition = """equal_to("VARIANTE", "RAPIDE")""",
                                                      SYME  =SIMP(statut='f', typ='TXM', defaut="OUI",
                                                                  into=("OUI","NON"),),
                                   ),
                                ),
# le choix du type de contact implique aussi celui de frottement
                TYPE_CONT  = SIMP(statut='f', typ='TXM', defaut="UNILATERAL",
                                    into=("UNILATERAL", "BILATERAL","COLLE"),
                                    fr=tr("Choix d'un modèle de contact" ),),

# coefficient de nitche, pénalisation ou augmentation en fonction de la méthode de contact
                COEF_CONT   =SIMP(statut='f',typ='R'  ,defaut=100.0, val_min=0.0),

# Pairing options (for segment to segment contact)
                # "COLLOCATION" old robuste, "RAPIDE" new => à discuter avec Guillaume
                APPARIEMENT     =SIMP(statut='f',typ='TXM',defaut="MORTAR", into=("MORTAR")),
                DIST_APPA       =SIMP(statut='f',typ='R'  ,defaut=-1.0),

                GROUP_MA_MAIT   =SIMP(statut='o',typ=grma ,max=1),
                GROUP_MA_ESCL   =SIMP(statut='o',typ=grma ,max=1),

# Managing boundary conditions with contact (slave side) que pour la méthode LAGRANGE
                b_cl            = BLOC(condition = """is_in("ALGO_CONT", ('LAGRANGIEN',))""",
                                       SANS_GROUP_MA   =SIMP(statut='f',typ=grma ,validators=NoRepeat(),max='**'),
                                        ),
                CONTACT_INIT    =SIMP(statut='f',typ='TXM',defaut="INTERPENETRE",
                                        into=("OUI","INTERPENETRE","NON"),),
                SEUIL_INIT      = SIMP(statut='f',typ='R'),

# Add suppl. gaps
                DIST_SUPP       =SIMP(statut='f',typ=(fonction_sdaster,nappe_sdaster,formule)),
                # À vérifier si l'algo marche pour POUTRE
                DIST_POUTRE     =SIMP(statut='f',typ='TXM',defaut="NON", into=("OUI","NON")),
                DIST_COQUE      =SIMP(statut='f',typ='TXM',defaut="NON", into=("OUI","NON")),
                b_cara=BLOC(condition = """equal_to("DIST_POUTRE", 'OUI') or equal_to("DIST_COQUE", 'OUI')""",
                            CARA_ELEM     =SIMP(statut='o',typ=(cara_elem) ),
                            ),
                ), # fin ZONE
        ), # fin BLOC


# ----- définition mot-clé facteur avec frottement
        b_zone_fric     = BLOC(condition = """equal_to("FROTTEMENT", 'OUI') """,
                ZONE = FACT(statut='o', max='**',
# VERIFICATION DE L'ORIENTATION ET DE LA COHERENCE DES NORMALES
                VERI_NORM       =SIMP(statut='f', typ='TXM', defaut="OUI", into=("OUI","NON"),
                                 fr=tr("Vérification de l'orientation (sortante) des normales aux surfaces"),),
# Method for contact
                ALGO_CONT   = SIMP(statut='f',typ='TXM',defaut="LAGRANGIEN",
                                        into=("LAGRANGIEN","NITSCHE","PENALISATION",),),

                b_algo_cont = BLOC(condition = """is_in("ALGO_CONT", ('LAGRANGIEN', 'NITSCHE'))""",
                                   VARIANTE  =SIMP(statut='f', typ='TXM', defaut="ROBUSTE",
                                                into=("RAPIDE","ROBUSTE"),
                                                fr=tr("Choix de la variante des formulations du contact" ),),
                                   b_vari_syme = BLOC(condition = """equal_to("VARIANTE", "RAPIDE")""",
                                                      SYME  =SIMP(statut='f', typ='TXM', defaut="OUI",
                                                                  into=("OUI","NON"),),
                                   ),
                                ),
# le choix du type de contact implique aussi celui de frottement
                TYPE_CONT  = SIMP(statut='f', typ='TXM', defaut="UNILATERAL",
                                    into=("UNILATERAL", "BILATERAL","COLLE"),
                                    fr=tr("Choix d'un modèle de contact" ),),

# coefficient de nitche, pénalisation ou augmentation en fonction de la méthode de contact
                COEF_CONT   =SIMP(statut='f',typ='R'  ,defaut=100.0, val_min=0.0),

# Pairing options (for segment to segment contact)
                # "COLLOCATION" old robuste, "RAPIDE" new => à discuter avec Guillaume
                APPARIEMENT     =SIMP(statut='f',typ='TXM',defaut="MORTAR", into=("MORTAR")),
                DIST_APPA       =SIMP(statut='f',typ='R'  ,defaut=-1.0),

                GROUP_MA_MAIT   =SIMP(statut='o',typ=grma ,max=1),
                GROUP_MA_ESCL   =SIMP(statut='o',typ=grma ,max=1),

# Managing boundary conditions with contact (slave side) que pour la méthode LAGRANGE
                b_cl            = BLOC(condition = """is_in("ALGO_CONT", ('LAGRANGIEN',))""",
                                       SANS_GROUP_MA   =SIMP(statut='f',typ=grma ,validators=NoRepeat(),max='**'),
                                        ),
                CONTACT_INIT    =SIMP(statut='f',typ='TXM',defaut="INTERPENETRE",
                                        into=("OUI","INTERPENETRE","NON"),),
                SEUIL_INIT      = SIMP(statut='f',typ='R'),

# Add suppl. gaps
                DIST_SUPP       =SIMP(statut='f',typ=(fonction_sdaster,nappe_sdaster,formule)),
                # À vérifier si l'algo marche pour POUTRE
                DIST_POUTRE     =SIMP(statut='f',typ='TXM',defaut="NON", into=("OUI","NON")),
                DIST_COQUE      =SIMP(statut='f',typ='TXM',defaut="NON", into=("OUI","NON")),
                b_cara=BLOC(condition = """equal_to("DIST_POUTRE", 'OUI') or equal_to("DIST_COQUE", 'OUI')""",
                            CARA_ELEM     =SIMP(statut='o',typ=(cara_elem) ),
                            ),

# FROTTEMENT

                COEF_FROT   =SIMP(statut='f',typ='R'  ,defaut=100.0, val_min=0.0),

                b_zone_lagr     = BLOC(condition = """equal_to("ALGO_CONT", 'LAGRANGIEN')""",
                        ALGO_FROT  =SIMP(statut='f',typ='TXM',defaut="LAGRANGIEN",
                                 into=("LAGRANGIEN",),),
                ),

                b_cont_nits = BLOC(condition = """equal_to("ALGO_CONT", 'NITSCHE')""",
                        ALGO_FROT       =SIMP(statut='f',typ='TXM',defaut="NITSCHE",
                                        into=("NITSCHE",),),
                ),

                b_cont_pena = BLOC(condition = """equal_to("ALGO_CONT", 'PENALISE')""",
                        ALGO_FROT       =SIMP(statut='f',typ='TXM',defaut="PENALISE",
                                        into=("PENALISE",),),
                ),

                b_from_frot_colle = BLOC(condition = """equal_to("TYPE_CONT", 'COLLE')""",
                            TYPE_FROT  =SIMP(statut='f', typ='TXM', defaut="COLLE", into=("COLLE"),
                               fr=tr("Choix d'un modèle de frottement" ),),
                                ),

                b_form_frot = BLOC(condition = """is_in("TYPE_CONT", ('UNILATERAL', 'BILATERAL'))""",
                            TYPE_FROT  =SIMP(statut='f', typ='TXM', defaut="SANS",
                                            into=("SANS", "TRESCA","COULOMB"),
                                            fr=tr("Choix d'un modèle de frottement" ),),
                            # les valeurs par défaut ???
                            b_tresca  = BLOC(condition = """equal_to("TYPE_FROT", "TRESCA")""",
                                        TRESCA  = SIMP(statut='o',typ='R'),),
                            b_coulomb = BLOC(condition = """equal_to("TYPE_FROT", "COULOMB")""",
                                        COULOMB  = SIMP(statut='o',typ='R'),),
                                ),

                ), # fin ZONE
        ), # fin BLOC

) #fin OPER


# Chantier  ADAPTATION à faire (2023)

# Ajouter une methode python pour activer cette option - doc dans la docstring
# -------------------- No contact solving (only pairing)  - Au niveau de la zone -----------------------------------------------------------------------------------
        # RESOLUTION       =SIMP(statut='f',typ='TXM',defaut="OUI",into=("OUI","NON")),
        #                         b_verif=BLOC(condition = """equal_to("RESOLUTION", 'NON') """,
        # TOLE_INTERP   = SIMP(statut='f',typ='R',defaut = 0.),),

#        ARRET DU CALCUL POUR LE MODE SANS RESOLUTION DU CONTACT - Option dans le python
        # STOP_INTERP   = SIMP(statut='f', typ='TXM', defaut="NON", into=("OUI","NON"),
        #                         fr=tr("Arrête le calcul dès qu'une interpénétration est détectée en mode RESOLUTION='NON'"),),




# Dans MECA_NON_LINE -> CONTACT

#         CONT_STAT_ELAS = SIMP(statut='f', typ='I',val_min=0, defaut = 0),

# PARAMETRE GENERAL : BOUCLE DE GEOMETRIE -> A mettre dans l'MKNL
#         ALGO_RESO_GEOM = SIMP(statut='f', typ='TXM',
#                                 into=("POINT_FIXE","NEWTON",), defaut="POINT_FIXE"),
#         b_algo_reso_geomNE = BLOC(condition = """equal_to("ALGO_RESO_GEOM", 'NEWTON')""",
#                                     RESI_GEOM = SIMP(statut='f',typ='R',defaut=1e-6),
#                                   ),
#         b_algo_reso_geomPF = BLOC(condition = """equal_to("ALGO_RESO_GEOM", 'POINT_FIXE')""",
#                                     REAC_GEOM = SIMP(statut='f',
#                                                    typ='TXM',
#                                                    into=("AUTOMATIQUE","CONTROLE","SANS",),
#                                                    defaut="AUTOMATIQUE"),
#                                     b_automatique = BLOC(condition = """equal_to("REAC_GEOM", 'AUTOMATIQUE') """,
#                                        ITER_GEOM_MAXI = SIMP(statut='f',typ='I',defaut=10),
#                                        RESI_GEOM      = SIMP(statut='f',typ='R',defaut=0.01)
#                                     ),
#                                     b_controle    = BLOC(condition = """equal_to("REAC_GEOM", 'CONTROLE') """,
#                                        NB_ITER_GEOM   = SIMP(statut='f',typ='I',defaut = 2),
#                                     ),
#            ),

# # PARAMETRE GENERAL : BOUCLE DE CONTACT (supprimé - Point fixe)
#   à vérifier si RESI_CONT a un sens pour la méthode NEWTON

# *****                    RESI_CONT=SIMP(statut='f',typ='R',defaut=-1,),



# # PARAMETRE GENERAL : BOUCLE DE FROTTEMENT
#         b_bouc_frot = BLOC(condition = """equal_to("FROTTEMENT", 'OUI') """,

#                                 b_algo_frot_geomNE = BLOC(condition = """equal_to("ALGO_RESO_GEOM", 'NEWTON')""",
#                                     ALGO_RESO_FROT = SIMP(statut='f', typ='TXM',
#                                                         into=("NEWTON",),
#                                                         defaut="NEWTON"),
#                                 ),
#                                 b_algo_frot_geomPF = BLOC(condition = """equal_to("ALGO_RESO_GEOM", 'POINT_FIXE')""",
#                                     ALGO_RESO_FROT = SIMP(statut='f', typ='TXM',
#                                                         into=("POINT_FIXE","NEWTON",),
#                                                         defaut="POINT_FIXE"),

#                                     b_algo_reso_frotPF = BLOC(condition = """equal_to("ALGO_RESO_FROT", 'POINT_FIXE')""",
#                                         ITER_FROT_MAXI = SIMP(statut='f',typ='I',defaut=10, val_min=0),
#                                     ),
#                                 ),

#                                 RESI_FROT      = SIMP(statut='f',typ='R',defaut=1.e-4),
#            ),

# AUTRES POINTES :
# ADAPTATION chantier 2023 => peut-être dans DEFI_LIST_INST
