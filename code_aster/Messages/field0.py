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

from ..Utilities import _

cata_msg = {
    10: _(
        """
Le maillage sur lequel s'appuie le modèle et le maillage du champ des contraintes ne sont pas les mêmes.
"""
    ),
    11: _(
        """Maille: %(k1)-8s - Nombre de points précédent: %(i1)d - Nombre de points courant: %(i2)d"""
    ),
    12: _(
        """Maille: %(k1)-8s - Nombre de composantes précédent: %(i1)d - Nombre de composantes courant: %(i2)d"""
    ),
    13: _(
        """
Erreur lors de la vérification de la cohérence entre les champs de contrainte.

On a affiché ci-dessus la liste des mailles pour lesquelles il y a des différences (nombre de points et nombre de composantes)
"""
    ),
}
