# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
from ...Supervis.ExecuteCommand import UserMacro
from .post_t_q_ops import post_t_q_ops

POST_T_Q_CATA=MACRO(nom="POST_T_Q",
                    op=OPS('Contrib.post_t_q_ops.post_t_q_ops'),
                    sd_prod=table_sdaster,
                    fr=tr("Calcul des facteurs d'intensité de contraintes en 2D et en 3D par "
                         "extrapolation des sauts de déplacements sur les lèvres de la fissure"),
                    reentrant='n',
                    #UIinfo={"groupes":("Post-traitements","Rupture",)},

           regles=(UN_PARMI('FISSURE','FOND_FISS'),),

         MODELISATION  =SIMP(statut='o',typ='TXM',
                             into=("3D","AXIS","D_PLAN","C_PLAN"),position='global',
                             fr=tr("Modélisation cohérente avec celle utilisée pour le calcul des déplacements")),
         FOND_FISS     =SIMP(statut='f',typ=fond_fissure),
         FISSURE       =SIMP(statut='f',typ=fiss_xfem),
         RESULTAT      =SIMP(statut='o',typ=(evol_elas,evol_noli,mode_meca),position='global',
                             fr=tr("Déplacement des noeuds de la lèvre supérieure et inférieure")),
         NB_NOEUD_COUPE=SIMP(statut='f',typ='I',defaut=5,val_min = 3),

#        bloc correspondant a la donnee du fond de fissure pour les fissures maillees
         b_fond_fiss   =BLOC (condition="FOND_FISS is not None",

              b_no_mod =BLOC (condition="AsType(RESULTAT)!= mode_meca",

                        EVOL_THER   = SIMP(statut='f',typ=(evol_ther),fr=tr("Température sur le fond de fissure")),
                        b_ref_3D    = BLOC (condition="MODELISATION=='3D' ",
                                  TYPE_MAILLAGE = SIMP(statut='f',typ='TXM',into=("LIBRE","REGLE"),defaut="REGLE"),
                                  ),
                             ),
              b_mod    =BLOC (condition="AsType(RESULTAT)== mode_meca and MODELISATION=='3D'",
                                  TYPE_MAILLAGE = SIMP(statut='f',typ='TXM',into=("REGLE",),defaut="REGLE"),
                             ),

                         NOEUD         = SIMP(statut='f',typ=no  ,validators=NoRepeat(),max='**'),
                         GROUP_NO      = SIMP(statut='f',typ=grno,validators=NoRepeat(),max='**'),
                         SANS_NOEUD    = SIMP(statut='f',typ=no  ,validators=NoRepeat(),max='**'),
                         SANS_GROUP_NO = SIMP(statut='f',typ=grno,validators=NoRepeat(),max='**'),
                         TOUT          = SIMP(statut='f',typ='TXM',into=("OUI",) ),
                         ABSC_CURV_MAXI=SIMP(statut='f',typ='R',
                             fr=tr("Distance maximum à partir du fond de fissure à utiliser pour le calcul")),
                         ),

#        bloc correspondant a la donnee de la fissure pour les fissures X-FEM
         b_fissure     =BLOC (condition="FISSURE is not None",
                         NB_POINT_FOND = SIMP(statut='f',typ='I' ,),
                         NUME_FOND     = SIMP(statut='f',typ='I',defaut=1),
                         ABSC_CURV_MAXI=SIMP(statut='f',typ='R',
                             fr=tr("Distance maximum à partir du fond de fissure à utiliser pour le calcul")),
                         ),

         MATER         =SIMP(statut='o',typ=mater_sdaster,
                             fr=tr("Matériau homogène et isotrope cohérent avec celui utilisé pour le calcul des déplacements")),

#         MAILLAGE        =SIMP(statut='o',typ=maillage_sdaster),


         PREC_VIS_A_VIS=SIMP(statut='f',typ='R',defaut=0.1),

         b_mod_meca  =BLOC (condition="AsType(RESULTAT)== mode_meca ",
         TOUT_ORDRE    =SIMP(statut='f',typ='TXM',into=("OUI",) ),
         NUME_ORDRE    =SIMP(statut='f',typ='I',validators=NoRepeat(),max='**'),
         LIST_ORDRE    =SIMP(statut='f',typ=listis_sdaster),
         TOUT_MODE       =SIMP(statut='f',typ='TXM',into=("OUI",) ),
         NUME_MODE       =SIMP(statut='f',typ='I',validators=NoRepeat(),max='**'),
         LIST_MODE       =SIMP(statut='f',typ=listis_sdaster),
         FREQ          =SIMP(statut='f',typ='R',validators=NoRepeat(),max='**'),
         LIST_FREQ     =SIMP(statut='f',typ=listr8_sdaster),
             b_acce_reel     =BLOC(condition="(FREQ is not None) or (LIST_FREQ is not None)",
               CRITERE         =SIMP(statut='f',typ='TXM',defaut="RELATIF",into=("RELATIF","ABSOLU") ),
               b_prec_rela=BLOC(condition="(CRITERE=='RELATIF')",
                   PRECISION       =SIMP(statut='f',typ='R',defaut= 1.E-6,),),
               b_prec_abso=BLOC(condition="(CRITERE=='ABSOLU')",
                   PRECISION       =SIMP(statut='o',typ='R',),),
             ),
             ),
         b_no_mod_meca  =BLOC (condition="AsType(RESULTAT)!= mode_meca ",
         TOUT_ORDRE    =SIMP(statut='f',typ='TXM',into=("OUI",) ),
         NUME_ORDRE    =SIMP(statut='f',typ='I',validators=NoRepeat(),max='**'),
         LIST_ORDRE    =SIMP(statut='f',typ=listis_sdaster),
         INST          =SIMP(statut='f',typ='R',validators=NoRepeat(),max='**'),
         LIST_INST     =SIMP(statut='f',typ=listr8_sdaster),
             b_acce_reel     =BLOC(condition="(INST is not None)or(LIST_INST is not None)",
               CRITERE         =SIMP(statut='f',typ='TXM',defaut="RELATIF",into=("RELATIF","ABSOLU") ),
               b_prec_rela=BLOC(condition="(CRITERE=='RELATIF')",
                   PRECISION       =SIMP(statut='f',typ='R',defaut= 1.E-6,),),
               b_prec_abso=BLOC(condition="(CRITERE=='ABSOLU')",
                   PRECISION       =SIMP(statut='o',typ='R',),),
             ),
             ),
         INFO          =SIMP(statut='f',typ='I',defaut=1,into=(1,2)),
         TITRE         =SIMP(statut='f',typ='TXM',max='**'),
)  ;

POST_T_Q = UserMacro("POST_T_Q", POST_T_Q_CATA,
                     post_t_q_ops)
