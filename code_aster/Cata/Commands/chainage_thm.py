# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2020 - EDF R&D - www.code-aster.org
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

# person_in_charge: sylvie.granet at edf.fr

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *


def chainage_thm_prod(self,TYPE_CHAINAGE,TYPE_RESU = None,**args) :
  if args.get('__all__'):
      return ([None, evol_varc, cham_no_sdaster],
              [None, corresp_2_mailla],
              [None, corresp_2_mailla],
              [None, corresp_2_mailla])

  if TYPE_CHAINAGE == "MECA_HYDR" : return evol_varc

  if TYPE_CHAINAGE == "HYDR_MECA" :
    if TYPE_RESU == "CHAM_NO" :
      return cham_no_sdaster
    elif TYPE_RESU == "EVOL_VARC" :
      return evol_varc

  if TYPE_CHAINAGE == "INIT" :
    matr_mh  = args['MATR_MH']
    matr_hm1 = args['MATR_HM1']
    matr_hm2 = args['MATR_HM2']

    self.type_sdprod(matr_mh,corresp_2_mailla)
    self.type_sdprod(matr_hm1,corresp_2_mailla)
    self.type_sdprod(matr_hm2,corresp_2_mailla)
    return None

  raise AsException("type de chainage THM non prevu")

CHAINAGE_THM=MACRO(nom="CHAINAGE_THM",
                   op=OPS('code_aster.MacroCommands.chainage_thm_ops.chainage_thm_ops'),
                   sd_prod=chainage_thm_prod,
                   reentrant='n',
                   docu="Ux.xx.xx",
                   fr=tr("Calcul des variables de commande pour le cha??nage THM"),

         TYPE_CHAINAGE  = SIMP(statut='o',typ='TXM',
                               into=("HYDR_MECA","MECA_HYDR","INIT",),
                               fr=tr("Sens du cha??nage ou initialisation des matrices de projection")),

         # Cas HYDR_MECA :

         b_hydr_meca    = BLOC(condition = """equal_to("TYPE_CHAINAGE", 'HYDR_MECA')""",fr=tr("Cha??nage hydraulique vers m??canique"),

             RESU_HYDR       = SIMP(statut='o',typ=resultat_sdaster,fr=tr("R??sultat hydraulique ?? cha??ner") ),
             MODELE_MECA     = SIMP(statut='o',typ=modele_sdaster  ,fr=tr("Mod??le d'arriv??e m??canique")),
             TYPE_RESU       = SIMP(statut='f',typ='TXM',into=("EVOL_VARC","CHAM_NO"),defaut="EVOL_VARC", ),
             MATR_HM1        = SIMP(statut='o',typ=corresp_2_mailla,),
             MATR_HM2        = SIMP(statut='o',typ=corresp_2_mailla,),

             b_type_resu     = BLOC(condition = """equal_to("TYPE_RESU", 'EVOL_VARC')""",fr=tr("Instant obligatoire si TYPE_RESU=EVOL_VARC"),
                                   INST = SIMP(statut='o',typ='R',max=1),

           ),),

         # Cas MECA_HYDR :

         b_meca_hydr    = BLOC(condition = """equal_to("TYPE_CHAINAGE", 'MECA_HYDR')""",fr=tr("Cha??nage m??canique vers hydraulique"),

             RESU_MECA       = SIMP(statut='o',typ=resultat_sdaster,fr=tr("R??sultat m??canique ?? cha??ner") ),
             MODELE_HYDR     = SIMP(statut='o',typ=modele_sdaster  ,fr=tr("Mod??le d'arriv??e hydraulique")),

             MATR_MH         = SIMP(statut='o',typ=corresp_2_mailla,),
             INST            = SIMP(statut='o',typ='R',max=1),
           ),

         # Cas INIT :

         b_init    = BLOC(condition = """equal_to("TYPE_CHAINAGE", 'INIT')""",fr=tr("Calcul des matrices de projection"),

             MODELE_MECA     = SIMP(statut='o',typ=modele_sdaster  ,fr=tr("Mod??le m??canique")),
             MODELE_HYDR     = SIMP(statut='o',typ=modele_sdaster  ,fr=tr("Mod??le hydraulique")),

             MATR_MH         = SIMP(statut='o',typ=CO,),
             MATR_HM1        = SIMP(statut='o',typ=CO,),
             MATR_HM2        = SIMP(statut='o',typ=CO,),
           ),

         INFO     = SIMP(statut='f',typ='I',defaut=1,into=( 1, 2 ) ),

) ;
