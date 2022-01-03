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

from ..Utilities import _

cata_msg = {

    4 : _("""Un des éléments esclave n'est pas du bon type.
 Pour le calcul de la normale, il faut que les éléments soient de la bonne dimension: des segments en 2D ou des faces en 3D."""),

    9 : _("""Il est interdit d'avoir deux mailles de type POI1 simultanément sur les deux surfaces en vis-à-vis."""),

   48 : _("""Il n'y a aucun noeud esclave à lier pour l'occurrence %(i1)d du mot clé LIAISON_MAIL.
   Peut-être que tous les noeuds esclaves ont déjà été éliminés dans des occurrences précédentes."""),

   49 : _(""" Pour le calcul de la normale sur le côté esclave, il faut donner des éléments de facette."""),

}
