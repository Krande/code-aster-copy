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

from ..Utilities import _

cata_msg = {

    1 : _("""Erreur lors de l'orientation des mailles de coque solide: il faut que le groupe de mailles volumiques ne contienne que des HEXA8 et que le groupe surfacique ne contienne que des QUAD4."""),

    2 : _("""Erreur lors de la création des mailles de coque solide: on doit orienter les hexaèdres, il faut utiliser GROUP_MA_SURF."""),

    3 : _("""Erreur lors de la création des mailles de coque solide: il n'y a que des pentaèdres, il est inutile d'utiliser GROUP_MA_SURF."""),
}
