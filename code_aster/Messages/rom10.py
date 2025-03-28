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

# Messages for mode management in ROM

from ..Utilities import _

cata_msg = {
    1: _("""La base contient des champs de type %(k1)s."""),
    2: _("""Un mode contient %(i1)d équations."""),
    10: _("""On ne peut définir des modes que sur des maillages tridimensionnels."""),
    11: _("""Le modèle doit être le même sur tous les modes."""),
    12: _("""On ne peut définir des modes que sur des modèles tridimensionnels."""),
}
