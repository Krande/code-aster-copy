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

EXTR_COUPE = MACRO(
    nom="EXTR_COUPE",
    op=OPS("code_aster.MacroCommands.extr_coupe_ops.extr_coupe_ops"),
    sd_prod=cham_gd_sdaster,
    fr=tr("Extraction de résultats sur un ensemble de coupes"),
    RESULTAT=SIMP(
        statut="o",
        typ=(evol_elas, evol_noli, evol_ther, evol_elas_dict, evol_ther_dict, evol_noli_dict),
    ),
    COUPES=SIMP(statut="o", typ=table_sdaster),
    NOM_CHAM=SIMP(statut="f", typ="TXM", validators=NoRepeat(), max="**", into=C_NOM_CHAM_INTO()),
    REPERE=SIMP(statut="f", typ="TXM", defaut="GLOBAL", into=("GLOBAL", "LOCAL")),
    LINEARISATION=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
    b_linearisation=BLOC(
        condition="""equal_to("LINEARISATION", "NON")""",
        REPARTITION=SIMP(
            statut="f", typ="TXM", defaut="UNIFORME", into=("UNIFORME", "GAUSSIENNE", "UTILISATEUR")
        ),
        b_repartition_utilisateur=BLOC(
            condition="""equal_to("REPARTITION", 'UTILISATEUR')""",
            FONCTION=SIMP(statut="o", typ=formule),
        ),  # fin bloc_repartition_utilisateur
        b_repartition_gaussien=BLOC(
            condition="""equal_to("REPARTITION", 'GAUSSIENNE')""",
            PARAM=SIMP(statut="o", typ="R", min=2, max=2),
        ),  # fin bloc_repartition_gaussien
    ),  # fin bloc_linearisation
)  # fin MACRO
