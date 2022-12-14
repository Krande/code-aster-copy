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


# definition of macro return types (both in-out arguments and return types)
def combinaison_ferraillage_prod ( self, **args ):
    return mult_elas

# definition of macro catalogue
COMBINAISON_FERRAILLAGE = MACRO(
    nom = 'COMBINAISON_FERRAILLAGE',
    op = OPS('code_aster.MacroCommands.combinaison_ferraillage_ops.combinaison_ferraillage_ops'),
    fr = tr("COMBINAISON_FERRAILLAGE"),
    sd_prod = combinaison_ferraillage_prod,
    reentrant = 'o:RESULTAT',
    reuse = SIMP ( statut = 'c', typ = CO ),
    RESULTAT = SIMP ( statut = 'o', typ = mult_elas ),

    COMBINAISON   = FACT ( statut = 'o', min = 1 , max = '**',
        TYPE = SIMP ( statut = 'o', typ = 'TXM',
                      into = (
                               'ELS_QUASIPERMANENT',
                               'ELS_CARACTERISTIQUE',
                               'ELU_FONDAMENTAL',
                               'ELU_ACCIDENTEL',
                              )
                    ),
        regles = ( UN_PARMI ( 'NOM_CAS', 'NUME_ORDRE' ), ),
           NOM_CAS = SIMP ( statut = 'f', typ = 'TXM', min = 1, max = '**' ),
        NUME_ORDRE = SIMP ( statut = 'f', typ = 'I', min = 1, max = '**' ),
    ),

    CODIFICATION  = SIMP ( statut = 'f', typ = 'TXM', defaut = 'EC2' ,
                           into = ('EC2','BAEL91') ),

         b_EC2 = BLOC(condition = """ equal_to("CODIFICATION", 'EC2')""",
                          fr=tr("utilisation de l'eurocode 2"),
#          mot cl?? facteur r??p??table pour assigner les caract??ristiques locales par zones topologiques (GROUP_MA)
           AFFE               =FACT(statut='o',max='**',
             regles           =(UN_PARMI('TOUT','GROUP_MA','MAILLE'),),
             TOUT             =SIMP(statut='f',typ='TXM',into=("OUI",) ),
             GROUP_MA         =SIMP(statut='f',typ=grma,validators=NoRepeat(),max='**'),
             MAILLE           =SIMP(statut='c',typ=ma,validators=NoRepeat(),max='**'),
             TYPE_STRUCTURE   =SIMP(statut = 'o', typ = 'TXM',
                                    into = ('1D', '2D',),
                                    fr=tr("Type de Structure 1D ou 2D")),
             FERR_SYME        =SIMP(statut = 'f', typ = 'TXM', defaut='NON',
                                    into = ('OUI', 'NON',),
                                    fr=tr("Ferraillage sym??trique?")),
             SEUIL_SYME       =SIMP(statut = 'f', typ = 'R',
                                    fr=tr("Seuil de tol??rance pour le calcul d'un ferraillage sym??trique")),
             FERR_COMP        =SIMP(statut = 'f', typ = 'TXM', defaut='NON',
                                    into = ('OUI', 'NON',),
                                    fr=tr("Ferraillage de compression possible?")),
             EPURE_CISA       =SIMP(statut = 'f', typ = 'TXM', defaut='NON',
                                    into = ('OUI', 'NON',),
                                    fr=tr("Prise en compte de l'effort de traction suppl??mentaire du ?? l'effort tranchant et ?? la torsion?")),
             UNITE_CONTRAINTE =SIMP(statut='o',typ='TXM', into=("MPa","Pa"),
                                    fr=tr("Unit?? des contraintes du probl??me")),
             UNITE_DIMENSION  =SIMP(statut='o',typ='TXM', into=("mm","m"),
                                    fr=tr("Unit?? des dimensions du probl??me")),
             FERR_MIN         =SIMP(statut='f',typ='TXM', defaut='NON',
                                    into=('NON','OUI','CODE',),
                                    fr=tr("Prise en compte d'un ferraillage minimal?")),
             RHO_LONGI_MIN    =SIMP(statut='f',typ='R',
                                    fr=tr("Ratio de ferraillage longitudinal minimal en %")),
             RHO_TRNSV_MIN    =SIMP(statut='f',typ='R',
                                    fr=tr("Ratio de ferraillage transversal minimal en %")),
             c_2D             =BLOC(condition = """ equal_to("TYPE_STRUCTURE", '2D')""",
                                    fr=tr("d??finition des enrobages de la section 2D"),
                                    C_INF=SIMP(statut='o',typ='R',
                                    fr=tr("Enrobage des armatures inf??rieures pour la section 2D")),
                                    C_SUP=SIMP(statut='o',typ='R',
                                    fr=tr("Enrobage des armatures sup??rieures pour la section 2D")),),
             c_1D             =BLOC(condition = """ equal_to("TYPE_STRUCTURE", '1D')""",
                                    fr=tr("d??finition des enrobages de la section 1D"),
                                    C_INF_Y=SIMP(statut='o',typ='R',
                                    fr=tr("Enrobage des armatures inf??rieures suivant l'axe Y de la section 1D")),
                                    C_SUP_Y=SIMP(statut='o',typ='R',
                                    fr=tr("Enrobage des armatures sup??rieures suivant l'axe Y de la section 1D")),
                                    C_INF_Z=SIMP(statut='o',typ='R',
                                    fr=tr("Enrobage des armatures inf??rieures suivant l'axe Z de la section 1D")),
                                    C_SUP_Z=SIMP(statut='o',typ='R',
                                    fr=tr("Enrobage des armatures sup??rieures suivant l'axe Z de la section 1D")),),
             ALPHA_E          =SIMP(statut='f',typ='R',
                                    fr=tr("Coefficient d'??quivalence acier/b??ton (ELS, ELS_QP)")),
             RHO_ACIER        =SIMP(statut='f',typ='R', defaut=-1,
                                    fr=tr("Densit?? volumique des aciers")),
             FYK              =SIMP(statut='f',typ='R',
                                    fr=tr("Limite d'??lasticit?? caract??ristique dans l'acier")),
             EYS              =SIMP(statut='f',typ='R',
                                    fr=tr("Module d'Young de l'acier")),
             TYPE_DIAGRAMME   =SIMP(statut='f',typ = 'TXM',defaut='B2',
                                    into = ('B1', 'B2',),
                                    fr=tr("Type du diagramme Contrainte-Deformation ?? utiliser: B1 (Inclin??) ou B2 (Horizontal)")),
             GAMMA_S_FOND     =SIMP(statut='f',typ='R',
                                    fr=tr("Coefficient de s??curit?? sur la r??sistance de calcul des aciers pour la combinaison fondamentale")),
             GAMMA_S_ACCI     =SIMP(statut='f',typ='R',
                                    fr=tr("Coefficient de s??curit?? sur la r??sistance de calcul des aciers pour la combinaison accidentelle")),
             FCK              =SIMP(statut='f',typ='R',
                                    fr=tr("R??sistance caract??ristique du b??ton en compression ?? 28 jours")),
             GAMMA_C_FOND     =SIMP(statut='f',typ='R',
                                    fr=tr("Coefficient de s??curit?? sur la r??sistance de calcul du b??ton pour la combinaison fondamentale")),
             GAMMA_C_ACCI     =SIMP(statut='f',typ='R',
                                     fr=tr("Coefficient de s??curit?? sur la r??sistance de calcul du b??ton pour la combinaison fondamentale")),           
             ALPHA_CC         =SIMP(statut='f',typ='R',defaut=1.,
                                    fr=tr("Coefficient de s??curit?? sur la r??sistance de calcul du b??ton en compression (ELU)")),
             SIGS_ELS         =SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement des aciers (ELS)")),
             sigc_2D          =BLOC(condition = """ equal_to("TYPE_STRUCTURE", '2D')""",
                                    fr=tr("d??finition des contraintes ultimes de dimensionnement du b??ton de la section 2D"),
                                    SIGC_INF_ELS=SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement du b??ton en fibre inf??rieure pour la section 2D (ELS)")),
                                    SIGC_SUP_ELS=SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement du b??ton en fibre sup??rieure pour la section 2D (ELS)")),),
             sigc_1D          =BLOC(condition = """ equal_to("TYPE_STRUCTURE", '1D')""",
                                    fr=tr("d??finition des contraintes ultimes de dimensionnement du b??ton de la section 1D"),
                                    SIGC_INF_Y_ELS=SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement du b??ton en fibre inf??rieure suivant l'axe Y de la section 1D (ELS)")),
                                    SIGC_SUP_Y_ELS=SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement du b??ton en fibre sup??rieure suivant l'axe Y de la section 1D (ELS)")),
                                    SIGC_INF_Z_ELS=SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement du b??ton en fibre inf??rieure suivant l'axe Z de la section 1D (ELS)")),
                                    SIGC_SUP_Z_ELS=SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement du b??ton en fibre sup??rieure suivant l'axe Z de la section 1D (ELS)")),),
             wmax_2D          =BLOC(condition = """ equal_to("TYPE_STRUCTURE", '2D')""",
                                    fr=tr("d??finition des ouvertures des fissures maximales de la section 2D"),
                                    WMAX_INF=SIMP(statut='f',typ='R',
                                    fr=tr("Ouverture maximale des fissures en face inf??rieure de la section 2D (ELS_QP)")),
                                    WMAX_SUP=SIMP(statut='f',typ='R',
                                    fr=tr("Ouverture maximale des fissures en face sup??rieure de la section 2D (ELS_QP)")),),
             wmax_1D          =BLOC(condition = """ equal_to("TYPE_STRUCTURE", '1D')""",
                                    fr=tr("d??finition des ouvertures des fissures maximales de la section 1D"),                                
                                    WMAX_INF_Y=SIMP(statut='f',typ='R',
                                    fr=tr("Ouverture maximale des fissures en face inf??rieure suivant l'axe Y de la section 1D (ELS_QP)")),
                                    WMAX_SUP_Y=SIMP(statut='f',typ='R',
                                    fr=tr("Ouverture maximale des fissures en face sup??rieure suivant l'axe Y de la section 1D (ELS_QP)")),
                                    WMAX_INF_Z=SIMP(statut='f',typ='R',
                                    fr=tr("Ouverture maximale des fissures en face inf??rieure suivant l'axe Z de la section 1D (ELS_QP)")),
                                    WMAX_SUP_Z=SIMP(statut='f',typ='R',
                                    fr=tr("Ouverture maximale des fissures en face sup??rieure suivant l'axe Z de la section 1D (ELS_QP)")),),
             SIGC_ELS_QP      =SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement du b??ton (ELS_QP)")),
             KT               =SIMP(statut='f',typ='R',
                                    fr=tr("Coefficient de dur??e de chargement (ELS_QP)")),
             PHI_INF_X        =SIMP(statut='f',typ='R',
                                    fr=tr("Diam??tre approximatif des armatures inf??rieures suivant l'axe X (ELS_QP)")),
             PHI_SUP_X        =SIMP(statut='f',typ='R',
                                    fr=tr("Diam??tre approximatif des armatures sup??rieures suivant l'axe X (ELS_QP)")),
             PHI_INF_Y        =SIMP(statut='f',typ='R',
                                    fr=tr("Diam??tre approximatif des armatures inf??rieures suivant l'axe Y (ELS_QP)")),
             PHI_SUP_Y        =SIMP(statut='f',typ='R',
                                    fr=tr("Diam??tre approximatif des armatures sup??rieures suivant l'axe Y (ELS_QP)")),
             PHI_INF_Z        =SIMP(statut='f',typ='R',
                                    fr=tr("Diam??tre approximatif des armatures inf??rieures suivant l'axe Z (ELS_QP)")),
             PHI_SUP_Z        =SIMP(statut='f',typ='R',
                                    fr=tr("Diam??tre approximatif des armatures sup??rieures suivant l'axe Z (ELS_QP)")),
             UTIL_COMPR       =SIMP(statut='f',typ='TXM',defaut='NON', into=("OUI","NON"),
                                    fr=tr("Prise en compte de la compression pour les aciers transversaux")),
             CLASSE_ACIER     =SIMP(statut='f',typ='TXM',defaut='B', into=("A","B","C"),
                                    fr=tr("Classe de ductilit?? des aciers")),
             b_iconst_ec2 = BLOC(condition = """ greater_than("RHO_ACIER", 0)""",
                                 fr=tr("Calcul du crit??re de difficult?? de b??tonnage si RHO_ACIER > 0"),
                ALPHA_REINF      =SIMP(statut='f',typ='R',defaut=1,
                                        fr=tr("Coefficient de pond??ration du ratio de densit?? d'acier par m??tre cube de b??ton")),
                ALPHA_SHEAR      =SIMP(statut='f',typ='R',defaut=1,
                                        fr=tr("Coefficient de pond??ration du ratio de densit?? d'acier d'effort tranchant")),
                ALPHA_STIRRUPS   =SIMP(statut='f',typ='R',defaut=1,
                                        fr=tr("Coefficient de pond??ration du ratio de longueur des ??pingles d'acier effort tranchant")),
                RHO_CRIT         =SIMP(statut='f',typ='R',defaut=150,
                                        fr=tr("Densit?? volumique d'armature critique")),
                DNSTRA_CRIT      =SIMP(statut='f',typ='R',defaut=0.006,
                                        fr=tr("Ferraillage d'effort tranchant critique")),
                L_CRIT           =SIMP(statut='f',typ='R',defaut=1,
                                        fr=tr("Longueur critique des epingle d'aciers d'effort tranchant")),
                ),
             ),),
                
         b_BAEL91 = BLOC(condition = """ equal_to("CODIFICATION", 'BAEL91')""",
                          fr=tr("utilisation du BAEL91"),
#          mot cl?? facteur r??p??table pour assigner les caract??ristiques locales par zones topologiques (GROUP_MA)
           AFFE               =FACT(statut='o',max='**',
             regles           =(UN_PARMI('TOUT','GROUP_MA','MAILLE'),),
             TOUT             =SIMP(statut='f',typ='TXM',into=("OUI",) ),
             GROUP_MA         =SIMP(statut='f',typ=grma,validators=NoRepeat(),max='**'),
             MAILLE           =SIMP(statut='c',typ=ma,validators=NoRepeat(),max='**'),
             TYPE_STRUCTURE   =SIMP(statut = 'o', typ = 'TXM',
                                    into = ('1D', '2D',),
                                    fr=tr("Type de Structure 1D ou 2D")),
             FERR_SYME        =SIMP(statut = 'f', typ = 'TXM', defaut='NON',
                                    into = ('OUI', 'NON',),
                                    fr=tr("Ferraillage sym??trique?")),
             SEUIL_SYME       =SIMP(statut = 'f', typ = 'R',
                                    fr=tr("Seuil de tol??rance pour le calcul d'un ferraillage sym??trique")),
             FERR_COMP        =SIMP(statut = 'f', typ = 'TXM', defaut='NON',
                                    into = ('OUI', 'NON',),
                                    fr=tr("Ferraillage de compression possible?")),
             EPURE_CISA       =SIMP(statut = 'f', typ = 'TXM', defaut='NON',
                                    into = ('OUI', 'NON',),
                                    fr=tr("Prise en compte de l'effort de traction suppl??mentaire du ?? l'effort tranchant et ?? la torsion?")),
             UNITE_CONTRAINTE =SIMP(statut='o',typ='TXM', into=("MPa","Pa"),
                                    fr=tr("Unit?? des contraintes du probl??me")),
             UNITE_DIMENSION  =SIMP(statut='o',typ='TXM', into=("mm","m"),
                                    fr=tr("Unit?? des dimensions du probl??me")),
             FERR_MIN         =SIMP(statut='f',typ='TXM', defaut='NON',
                                    into=('NON','OUI','CODE',),
                                    fr=tr("Prise en compte d'un ferraillage minimal?")),
             RHO_LONGI_MIN    =SIMP(statut='f',typ='R',
                                    fr=tr("Ratio de ferraillage longitudinal minimal en %")),
             RHO_TRNSV_MIN    =SIMP(statut='f',typ='R',
                                    fr=tr("Ratio de ferraillage transversal minimal en %")),
             c_2D             =BLOC(condition = """ equal_to("TYPE_STRUCTURE", '2D')""",
                                    fr=tr("d??finition des enrobages de la section 2D"),
                                    C_INF=SIMP(statut='o',typ='R',
                                    fr=tr("Enrobage des armatures inf??rieures pour la section 2D")),
                                    C_SUP=SIMP(statut='o',typ='R',
                                    fr=tr("Enrobage des armatures sup??rieures pour la section 2D")),),
             c_1D             =BLOC(condition = """ equal_to("TYPE_STRUCTURE", '1D')""",
                                    fr=tr("d??finition des enrobages de la section 1D"),
                                    C_INF_Y=SIMP(statut='o',typ='R',
                                    fr=tr("Enrobage des armatures inf??rieures suivant l'axe Y de la section 1D")),
                                    C_SUP_Y=SIMP(statut='o',typ='R',
                                    fr=tr("Enrobage des armatures sup??rieures suivant l'axe Y de la section 1D")),
                                    C_INF_Z=SIMP(statut='o',typ='R',
                                    fr=tr("Enrobage des armatures inf??rieures suivant l'axe Z de la section 1D")),
                                    C_SUP_Z=SIMP(statut='o',typ='R',
                                    fr=tr("Enrobage des armatures sup??rieures suivant l'axe Z de la section 1D")),),
             N                =SIMP(statut='f',typ='R',
                                    fr=tr("Coefficient d'??quivalence acier/b??ton (ELS,ELS_QP)")),
             RHO_ACIER        =SIMP(statut='f',typ='R', defaut=-1,
                                   fr=tr("Densit?? volumique des aciers")),
             FE               =SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte admissible dans l'acier")),
             EYS              =SIMP(statut='f',typ='R',
                                    fr=tr("Module d'Young de l'acier")),
             TYPE_DIAGRAMME   =SIMP(statut='f',typ = 'TXM',
                                    into = ('B1', 'B2',),
                                    fr=tr("Type du diagramme Contrainte-Deformation ?? utiliser: B1 (Inclin??) ou B2 (Horizontal)")),
             GAMMA_S          =SIMP(statut='f',typ='R',
                                    fr=tr("Coefficient de s??curit?? sur la r??sistance de calcul des aciers ?? l'ELU")),
             FCJ              =SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte admissible dans le b??ton")),
             GAMMA_C          =SIMP(statut='f',typ='R',
                                    fr=tr("Coefficient de s??curit?? sur la r??sistance de calcul du b??ton ?? l'ELU")),
             ALPHA_CC         =SIMP(statut='f',typ='R',defaut=0.85,
                                    fr=tr("Coefficient de s??curit?? sur la r??sistance de calcul du b??ton en compression (ELU)")),
             SIGS_ELS         =SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement des aciers ?? l'ELS")),
             sigc_2D          =BLOC(condition = """ equal_to("TYPE_STRUCTURE", '2D')""",
                                    fr=tr("d??finition des contraintes ultimes de dimensionnement du b??ton de la section 2D"),
                                    SIGC_INF_ELS=SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement du b??ton en fibre inf??rieure pour la section 2D (ELS)")),
                                    SIGC_SUP_ELS=SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement du b??ton en fibre sup??rieure pour la section 2D (ELS)")),),
             sigc_1D          =BLOC(condition = """ equal_to("TYPE_STRUCTURE", '1D')""",
                                    fr=tr("d??finition des contraintes ultimes de dimensionnement du b??ton de la section 1D"),
                                    SIGC_INF_Y_ELS=SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement du b??ton en fibre inf??rieure suivant l'axe Y de la section 1D (ELS)")),
                                    SIGC_SUP_Y_ELS=SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement du b??ton en fibre sup??rieure suivant l'axe Y de la section 1D (ELS)")),
                                    SIGC_INF_Z_ELS=SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement du b??ton en fibre inf??rieure suivant l'axe Z de la section 1D (ELS)")),
                                    SIGC_SUP_Z_ELS=SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement du b??ton en fibre sup??rieure suivant l'axe Z de la section 1D (ELS)")),),
             wmax_2D          =BLOC(condition = """ equal_to("TYPE_STRUCTURE", '2D')""",
                                    fr=tr("d??finition des ouvertures des fissures maximales de la section 2D"),
                                    WMAX_INF=SIMP(statut='f',typ='R',
                                    fr=tr("Ouverture maximale des fissures en face inf??rieure de la section 2D (ELS_QP)")),
                                    WMAX_SUP=SIMP(statut='f',typ='R',
                                    fr=tr("Ouverture maximale des fissures en face sup??rieure de la section 2D (ELS_QP)")),),
             wmax_1D          =BLOC(condition = """ equal_to("TYPE_STRUCTURE", '1D')""",
                                    fr=tr("d??finition des ouvertures des fissures maximales de la section 1D"),                                
                                    WMAX_INF_Y=SIMP(statut='f',typ='R',
                                    fr=tr("Ouverture maximale des fissures en face inf??rieure suivant l'axe Y de la section 1D (ELS_QP)")),
                                    WMAX_SUP_Y=SIMP(statut='f',typ='R',
                                    fr=tr("Ouverture maximale des fissures en face sup??rieure suivant l'axe Y de la section 1D (ELS_QP)")),
                                    WMAX_INF_Z=SIMP(statut='f',typ='R',
                                    fr=tr("Ouverture maximale des fissures en face inf??rieure suivant l'axe Z de la section 1D (ELS_QP)")),
                                    WMAX_SUP_Z=SIMP(statut='f',typ='R',
                                    fr=tr("Ouverture maximale des fissures en face sup??rieure suivant l'axe Z de la section 1D (ELS_QP)")),),
             SIGC_ELS_QP      =SIMP(statut='f',typ='R',
                                    fr=tr("Contrainte ultime de dimensionnement du b??ton (ELS_QP)")),
             KT               =SIMP(statut='f',typ='R',
                                    fr=tr("Coefficient de dur??e de chargement (ELS_QP)")),
             PHI_INF_X        =SIMP(statut='f',typ='R',
                                    fr=tr("Diam??tre approximatif des armatures inf??rieures suivant l'axe X (ELS_QP)")),
             PHI_SUP_X        =SIMP(statut='f',typ='R',
                                    fr=tr("Diam??tre approximatif des armatures sup??rieures suivant l'axe X (ELS_QP)")),
             PHI_INF_Y        =SIMP(statut='f',typ='R',
                                    fr=tr("Diam??tre approximatif des armatures inf??rieures suivant l'axe Y (ELS_QP)")),
             PHI_SUP_Y        =SIMP(statut='f',typ='R',
                                    fr=tr("Diam??tre approximatif des armatures sup??rieures suivant l'axe Y (ELS_QP)")),
             PHI_INF_Z        =SIMP(statut='f',typ='R',
                                    fr=tr("Diam??tre approximatif des armatures inf??rieures suivant l'axe Z (ELS_QP)")),
             PHI_SUP_Z        =SIMP(statut='f',typ='R',
                                    fr=tr("Diam??tre approximatif des armatures sup??rieures suivant l'axe Z (ELS_QP)")),
             b_iconst_bael = BLOC(condition = """ greater_than("RHO_ACIER", 0)""",
                                  fr=tr("Calcul du crit??re de difficult?? de b??tonnage si RHO_ACIER > 0"),
                ALPHA_REINF      =SIMP(statut='f',typ='R',defaut=1,
                                        fr=tr("Coefficient de pond??ration du ration de densit?? d'acier par m??tre cube de b??ton")),
                ALPHA_SHEAR      =SIMP(statut='f',typ='R',defaut=1,
                                        fr=tr("Coefficient de pond??ration du ration de densit?? d'acier d'effort tranchant")),
                ALPHA_STIRRUPS   =SIMP(statut='f',typ='R',defaut=1,
                                        fr=tr("Coefficient de pond??ration du ration de longueur des ??pingles d'acier d'effort tranchant")),
                RHO_CRIT         =SIMP(statut='f',typ='R',defaut=150,
                                        fr=tr("Densit?? volumique d'armature critique")),
                DNSTRA_CRIT      =SIMP(statut='f',typ='R',defaut=0.006,
                                        fr=tr("Ferraillage d'effort tranchant critique")),
                L_CRIT           =SIMP(statut='f',typ='R',defaut=1,
                                        fr=tr("Longueur critique des ??pingles d'aciers d'effort tranchant")),
                ),
             ),),
       )
