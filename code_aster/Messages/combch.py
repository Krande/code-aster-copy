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
 Le chargement %(k1)s n'est pas défini dans la table des coefficients
"""
    ),
    2: _(
        """
 La combinaison à la ligne %(i1)d dans la table applique plus de 2 chargements thermiques.
 Merci de vérifier votre table des combinaisons.
"""
    ),
    3: _(
        """
 S'il existe des chargements thermiques pour des éléments de structure MEMBRANE, 
 %(k1)s, %(k2)s, %(k3)s doivent être renseignés. 

 Car les chargements thermiques ne sont pas possible pour les éléments MEMBRANE.
 Il faut créer un modèle mécanique pour les calculs thermiques en remplaçant MEMBRANE par DKT.
"""
    ),
    4: _(
        """
  -----------------------------------------------------
 >>>>>>>>>> %(k1)s <<<<<<<<<<
 -----------------------------------------------------
"""
    ),
}
