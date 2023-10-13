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

from ..Utilities import _

cata_msg = {
    2: _(
        """
 Erreur dans ENV_CINE_YACS : problème avec VIS_A_VIS
"""
    ),
    8: _(
        """
 Erreur dans PROJ_CHAMP option COUPLAGE
     Nombre d'interfaces définies dans Code_Saturne : %(i1)i
     Nombre d'interfaces définies dans Code_Aster   : %(i2)i
     Vérifiez la cohérence entre :
       Code_Aster   : Le mot-clé GROUP_MA_IFS de la commande CALC_IFS_DNL
       Code_Saturne : La définition des structures dans USASTE.F
"""
    ),
    9: _(
        """
 Erreur dans PROJ_CHAMP option COUPLAGE : type de maille non reconnue : %(k1)s
"""
    ),
}
