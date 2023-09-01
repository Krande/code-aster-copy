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
from ...SD.sd_table import sd_table

CREA_COUPE = MACRO(
    nom="CREA_COUPE",
    op=OPS("code_aster.MacroCommands.crea_coupe_ops.crea_coupe_ops"),
    sd_prod=sd_table,
    fr=tr(
        """Création d'une table de coupes avec
                          des noeuds projetés sur la peau du maillage"""
    ),
    COUPE=SIMP(status="o", typ=table_sdaster),
    MAILLAGE=SIMP(status="o", typ=(maillage_sdaster, maillage_p)),
    # Nommage automatique des coupes
    NOM_AUTO=SIMP(status="o", typ="TXM", defaut="NON", into=("OUI", "NON")),
    b_nom_auto=BLOC(
        condition="""(equal_to("NOM_AUTO", 'OUI'))""",
        PREFIX=SIMP(
            status="o", typ="TXM", defaut="", fr=tr("préfixe à ajouter au début du nom des coupes")
        ),
        PREM_NUME=SIMP(
            status="f", typ="I", defaut=1, fr=tr("numéro à partir duquel la numérotation commence")
        ),
        PAS_NUME=SIMP(
            status="f",
            typ="I",
            defaut=1,
            fr=tr("incrément de numérotation entre deux coupes successives"),
        ),
    ),  # fin bloc b_nom_auto
)  # fin MACRO
