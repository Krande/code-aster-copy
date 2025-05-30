# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

from ..Language.DataStructure import UnitType
from ..Language.Syntax import FACT, SIMP


def C_AFFICHAGE():
    return FACT(
        statut="d",
        max=1,
        INFO_RESIDU=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
        INFO_TEMPS=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
        UNITE=SIMP(statut="f", typ=UnitType(), inout="out"),
        PAS=SIMP(statut="f", typ="I", val_min=1),
    )
