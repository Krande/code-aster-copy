# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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

# Pour la méthode d'adaptation du pas de temps

from ..Utilities import _

cata_msg = {
    1: _(
        """
  Adaptation du pas de temps.
"""
    ),
    2: _(
        """
    Pour la méthode d'adaptation de type <%(k1)s>, le pas de temps calculé vaut <%(r1)19.12e>.
"""
    ),
    3: _(
        """
    Pour la méthode d'adaptation de type <%(k1)s>, le critère n'est pas vérifié. Le pas de temps n'est pas adapté.
"""
    ),
    4: _(
        """
    Aucun critère d'adaptation n'est vérifié. On garde le pas de temps <%(r1)19.12e>.
"""
    ),
    5: _(
        """
    Sur tous les critères d'adaptation, le plus petit pas de temps vaut <%(r1)19.12e>.
"""
    ),
    6: _(
        """
    Après ajustement sur les points de passage obligatoires, le plus petit pas de temps vaut <%(r1)19.12e>.
"""
    ),
    10: _(
        """
    On maintient la découpe du pas de temps à <%(r1)19.12e>.
"""
    ),
    11: _(
        """
    La valeur du pas de temps retenu <%(r1)19.12e> est inférieure à PAS_MINI.
"""
    ),
    12: _(
        """
    La valeur du pas de temps <%(r1)19.12e> est supérieure à PAS_MAXI <%(r2)19.12e>.
    On limite le pas de temps à PAS_MAXI <%(r2)19.12e>.
"""
    ),
    13: _(
        """On a atteint le nombre maximal de pas de temps autorisé dans DEFI_LIST_INST. On arrête le calcul."""
    ),
}
