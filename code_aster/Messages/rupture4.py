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
    1: _(
        """
Résultat du test faux pour le produit vectoriel %(r1)f.
Il s'agit d'un test purement informatique. 
On continue le calcul. 
"""
    ),
    5: _(
        """
La valeur de NUME_ORDRE n'est pas valide.
"""
    ),
    6: _(
        """
La valeur de INST n'est pas valide.
"""
    ),
    7: _(
        """
Seule l'option GROUP_NO ou NB_POINT_FOND peut être sélectionnée.
"""
    ),
    8: _(
        """
Erreur sur le nombre de noeuds de réduction.
"""
    ),
    9: _(
        """
Le nombre de contours (NB_COUCHES) doit  être supérieur ou égal à 4.
"""
    ),
    10: _(
        """
Il est nécessaire de définir au moins une lèvre.
"""
    ),
    11: _(
        """
Le nombre de contours (NB_COUCHES) doit être inférieur à 20.
"""
    ),
    12: _(
        """
Le nombre de contours (NB_COUCHES) doit être inférieur à 10.
"""
    ),
    13: _(
        """
Le nombre de contours (NB_COUCHES) doit  être supérieur ou égal à 2.
"""
    ),
}
