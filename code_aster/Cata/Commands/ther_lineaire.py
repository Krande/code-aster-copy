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

# person_in_charge: mickael.abbas at edf.fr

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *

THER_LINEAIRE=MACRO(nom="THER_LINEAIRE",
                   op=OPS("code_aster.MacroCommands.ther_lineaire_ops.ther_lineaire_ops"),
                   sd_prod=evol_ther,
                   reentrant='f:RESULTAT',
                   fr=tr("Résoudre un problème thermique linéaire stationnaire ou transitoire"),
                   reuse=SIMP(statut='c', typ=CO),
                   RESULTAT        =SIMP(statut='f',typ=evol_ther,
                                         fr=tr("Objet qui sera enrichi des nouveaux instants calculés")),
                   MODELE          =SIMP(statut='o',typ=modele_sdaster),
                   CHAM_MATER      =SIMP(statut='o',typ=cham_mater),
                   CARA_ELEM       =SIMP(statut='f',typ=cara_elem),
                   EXCIT           =FACT(statut='o',max='**',
                                         CHARGE          =SIMP(statut='o',typ=(char_ther,char_cine_ther)),
                                         FONC_MULT       =SIMP(statut='f',typ=(fonction_sdaster,nappe_sdaster,formule)),
                   ),
                   ETAT_INIT       =FACT(statut='f', max=1,
                                         regles=(EXCLUS('STATIONNAIRE','EVOL_THER','CHAM_NO','VALE'),),
                                         STATIONNAIRE    =SIMP(statut='f',typ='TXM',into=("OUI",)),
                                         EVOL_THER       =SIMP(statut='f',typ=evol_ther),
                                         CHAM_NO         =SIMP(statut='f',typ=cham_no_sdaster),
                                         VALE            =SIMP(statut='f',typ='R'),
                                         NUME_ORDRE      =SIMP(statut='f',typ='I'),
                                         INST            =SIMP(statut='f',typ='R'),
                                         CRITERE         =SIMP(statut='f',typ='TXM',defaut="RELATIF",into=("RELATIF","ABSOLU") ),
                                         b_prec_rela=BLOC(condition="""(equal_to("CRITERE", 'RELATIF'))""",
                                                          PRECISION       =SIMP(statut='f',typ='R',defaut= 1.E-6,),),
                                         b_prec_abso=BLOC(condition="""(equal_to("CRITERE", 'ABSOLU'))""",
                                                          PRECISION       =SIMP(statut='o',typ='R',),),
                   ),
                   #-------------------------------------------------------------------
                   INCREMENT       =C_INCREMENT('THERMIQUE'),
                   #-------------------------------------------------------------------
                   #        Catalogue commun SOLVEUR
                   SOLVEUR         =C_SOLVEUR('THER_LINEAIRE'),
                   #-------------------------------------------------------------------
                   PARM_THETA      =SIMP(statut='f',typ='R',defaut= 0.57, val_min=0., val_max = 1.),
                   #-------------------------------------------------------------------
                   ARCHIVAGE       =C_ARCHIVAGE(),
                   #-------------------------------------------------------------------
                   TITRE           =SIMP(statut='f',typ='TXM'),
                   INFO            =SIMP(statut='f',typ='I',into=(1,2)),
                   translation={
                       "THER_LINEAIRE": "Linear thermal analysis",
                   }
)
