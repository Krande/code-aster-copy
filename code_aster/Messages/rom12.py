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

# Messages for base management in ROM

from ..Utilities import _

cata_msg = {
    1: _("""Construction de la matrice des modes de dimension [%(i1)d, %(i2)d]."""),
    2: _("""Sauvegarde de la base avec %(i1)d modes."""),
    3: _("""Création de la base pour les modes réduits avec %(i1)d numéros d'ordre."""),
    10: _("""Paramètres de la base:"""),
    11: _("""La base contient %(i1)d modes."""),
    12: _("""La base contient des modes linéiques."""),
    13: _("""Les modes linéiques ont pour axe: %(k1)s"""),
    14: _("""Les modes linéiques ont pour section de référence le GROUP_MA %(k1)s ."""),
    15: _("""La base contient des modes volumiques."""),
    16: _("""La base a été construite avec %(i1)d clichés."""),
    17: _("""La base contient les modes suivants:"""),
}
