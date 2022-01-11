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

MECA_NON_LINE = MACRO(nom="MECA_NON_LINE",
                      op=OPS(
                          "code_aster.MacroCommands.meca_non_line_ops.meca_non_line_ops"),
                      sd_prod=evol_noli,
                      fr=tr("Calcul de l'évolution mécanique ou thermo-hydro-mécanique couplée, en quasi-statique,"
                            " d'une structure en non linéaire"),
                      reentrant='f:RESULTAT',
                      reuse=SIMP(statut='c', typ=CO),
                      # -------------------------------------------------------------------
                      RESULTAT=SIMP(statut='f', typ=evol_noli,
                                    fr=tr("Objet qui sera enrichi des nouveaux instants calculés")),
                      # -------------------------------------------------------------------
                      MODELE=SIMP(statut='o', typ=modele_sdaster),
                      # -------------------------------------------------------------------
                      CHAM_MATER=SIMP(statut='o', typ=cham_mater),
                      # -------------------------------------------------------------------
                      EXCIT=FACT(statut='f', max='**',
                                 CHARGE=SIMP(statut='o', typ=(
                                     char_meca, char_cine_meca)),
                                 FONC_MULT=SIMP(statut='f', typ=(
                                     fonction_sdaster, nappe_sdaster, formule)),
                                 TYPE_CHARGE=SIMP(statut='f', typ='TXM', defaut="FIXE_CSTE",
                                                  into=("FIXE_CSTE",)),
                                 ),
                      # -------------------------------------------------------------------
                      SCHEMA_THM=C_SCHEMA_THM(),
                      # -------------------------------------------------------------------
                      COMPORTEMENT=C_COMPORTEMENT('STAT_NON_LINE'),
                      # -------------------------------------------------------------------
                      ETAT_INIT=C_ETAT_INIT('STAT_NON_LINE', 'f'),
                      # -------------------------------------------------------------------
                      INCREMENT=C_INCREMENT('MECANIQUE'),
                      # -------------------------------------------------------------------
                      METHODE=SIMP(statut='f', typ='TXM', defaut="NEWTON", into=(
                          "NEWTON", "NEWTON_KRYLOV")),
                      b_meth_newton=BLOC(condition="""equal_to("METHODE", 'NEWTON') or equal_to("METHODE", 'NEWTON_KRYLOV')""",
                                         NEWTON=C_NEWTON(),
                                         ),
                      # -------------------------------------------------------------------
                      CONVERGENCE=C_CONVERGENCE(),
                      # -------------------------------------------------------------------
                      SOLVEUR=C_SOLVEUR('STAT_NON_LINE'),
                      # -------------------------------------------------------------------
                      INFO=SIMP(statut='f', typ='I', into=(1, 2)),
                      )
