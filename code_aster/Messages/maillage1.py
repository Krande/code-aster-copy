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
    1: _(
        """
 Erreur de maillage :
   La maille %(k1)s de type %(k2)s est trop distordue.
   Le jacobien de la transformation géométrique n'a pas le même signe sur tous les
   points de Gauss.

 Risques & conseils :
   Le maillage a-t-il été produit par un mailleur ?
   La connectivité respecte-t-elle bien la convention Aster ?
"""
    ),
    2: _(
        """
Pour le noeud %(k1)s de la maille %(k2)s, la coordonnée X est négative (x=%(r1)G).
"""
    ),
    3: _(
        """
Pour une modélisation axisymétrique, la coordonnée X doit être positive, nulle ou
très faiblement négative ( > -1.d-6 * X_MAX)

 Conseils :
  * Vérifiez votre maillage.
  * Vous pouvez utiliser MODI_MAILLAGE / DEFORME pour repositionner votre maillage
    dans le demi espace  X >= 0
"""
    ),
    4: _(
        """
Il y a des noeuds doubles à l'interface de deux sous-domaines.
Il n'est pas possible de créer les joints avec.

 Conseils :
  * Supprimer les noeuds doubles.
  * Utiliser un nombre différents de processus MPI pour ne plus avoir de noeuds doubles
  à une interface.
"""
    ),
}
