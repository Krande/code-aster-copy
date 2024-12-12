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

# person_in_charge: francesco.bettonte at edf.fr

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *

LIRE_CORR_HOMO = MACRO(
    nom="LIRE_CORR_HOMO",
    op=OPS("code_aster.MacroCommands.MateHomo.lire_corr_homo_ops.lire_corr_ops"),
    sd_prod=table_sdaster,
    docu="UX.YZ.AB",
    reentrant="n",
    fr=tr("Lecture des correcteurs d'homogénéisation à partir d'un fichier MED."),
    UNITE=SIMP(statut="o", typ=UnitType(), inout="in"),
    CORR_MECA=SIMP(statut="f", typ=CO),
    CORR_THER=SIMP(statut="f", typ=CO),
)
