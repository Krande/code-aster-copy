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

#

from ..Utilities import _

cata_msg = {

     1 : _("""Type de maille non-supporté dans la projection. Il faut demander un développement."""),

    54 : _("""Il n'y a aucun noeud sur lesquels projeter."""),

    55 : _("""Il n'y a pas de mailles à projeter ou en correspondance.
 Dans le cas de l'utilisation de AFFE_CHAR_MECA / LIAISON_MAIL, les mailles maîtres
 doivent avoir la même dimension que l'espace de modélisation :
 - mailles volumiques pour un modèle 3D
 - mailles surfaciques pour un modèle 2D
"""),

    56 : _("""Le noeud %(k1)s n'a pas été trouvé lors de la projection. """),


}
