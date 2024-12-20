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

# For the command CALCUL and not for CALCUL.F90 !!!

from ..Utilities import _

cata_msg = {
    1: _(
        """
La table que l'on tente d'enrichir dans la commande CALCUL n'est pas du type attendu.
Elle n'a pas le bon nombre de paramètres.
"""
    ),
    2: _(
        """
La table que l'on tente d'enrichir dans la commande CALCUL n'est pas du type attendu.
Elle n'a pas les bons paramètres.
"""
    ),
    3: _(
        """
Le nom de la table donné par le mot-clef TABLE n'est pas le même que celui de la table produite par la commande CALCUL.
"""
    ),
    4: _(
        """
L'objet %(k1)s au numéro d'ordre %(i1)d existe déjà dans la table fournie.
On l'écrase pour le remplacer par le nouveau.
"""
    ),
    5: _(
        """
A cause des erreurs précédentes, le code s'arrête.
  Le champ des variables internes fourni à CALCUL n'est pas cohérent avec le comportement donné par le mot-clef COMPORTEMENT.
"""
    ),
    6: _(
        """
A cause des erreurs précédentes, le code s'arrête.
Le champ des contraintes fourni à CALCUL n'est pas cohérent à propos du nombre des sous-points.
"""
    ),
    7: _(
        """
Le modèle contient des éléments qui ne supportent pas les champs de type variables internes.
"""
    ),
    8: _(
        """
Le modèle contient des éléments qui ne supportent pas les champs de type variables internes.
On ne peut donc pas calculer les options les réclamant.
"""
    ),
    9: _(
        """
Le modèle contient des éléments qui ne supportent pas les champs de type variables internes.
Or les variables de commande ont des champs métallurgiques. On ne peut pas calculer la contribution de la métallurgie.
"""
    ),
    10: _(
        """
Le modèle contient des éléments XFEM.
La commande CALCUL ne sait pas traiter ce cas.
"""
    ),
    11: _(
        """
Le modèle contient des macro-éléments statiques.
La commande CALCUL ne sait pas traiter ce cas.
"""
    ),
    12: _(
        """
Pour calculer des options non-linéaires dans la commande CALCUL, il faut les déplacements, les contraintes et les variables internes précédents.
"""
    ),
    13: _(
        """
Pour calculer l'option FORC_NODA dans la commande CALCUL, il faut les déplacements et les contraintes précédents.
"""
    ),
    14: _(
        """
Le déplacement et l'incrément de déplacement ne sont pas définis sur la même numérotation. Soit le maillage est différent, soit ils ne reposent pas sur le même modèle, soit ils n'ont pas les mêmes conditions limites dualisées.
"""
    ),
}
