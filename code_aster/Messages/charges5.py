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
    2: _(
        """
 La charge de nom <%(k1)s> est en double.
"""
    ),
    3: _(
        """
 La charge de nom <%(k1)s> n'est pas autorisée dans la commande <%(k2)s>.
"""
    ),
    4: _(
        """
 La charge de nom <%(k1)s> se base sur un PHENOMENE non pris en charge dans la commande.
"""
    ),
    5: _(
        """
 La charge de nom <%(k1)s> se base sur un MODELE différent de la commande.
"""
    ),
    6: _(
        """
 Le mot-clé MODELE est nécessaire si CHARGE est présent dans une occurrence de EXCIT.
"""
    ),
    7: _(
        """
 La charge de nom <%(k1)s> ne peut pas être suiveuse.
 Les charges de type Dirichlet par élimination (AFFE_CHAR_CINE) ne peuvent pas être des charges suiveuses.
"""
    ),
    8: _(
        """
 La charge de nom <%(k1)s> ne peut pas être DIDI.
 Les charges de type Dirichlet par élimination (AFFE_CHAR_CINE) ne peuvent pas être des charges DIDI.
"""
    ),
    9: _(
        """
 La charge de nom <%(k1)s> ne peut pas être pilotée.
 Les charges de type Dirichlet par élimination (AFFE_CHAR_CINE) ne peuvent pas être des charges pilotées.
"""
    ),
    10: _(
        """
 La charge de nom <%(k1)s> ne peut pas être suiveuse.
 Les charges de type Dirichlet ne peuvent pas être des charges suiveuses.
"""
    ),
    11: _(
        """
 La charge de nom <%(k1)s> ne peut pas être pilotée.
 Les charges de type EVOL_CHAR ne peuvent pas être des charges pilotées.
"""
    ),
    12: _(
        """
 La charge de nom <%(k1)s> ne peut pas être associée à une fonction multiplicatrice complexe.
 Les charges de type Dirichlet par élimination (AFFE_CHAR_CINE) ne peuvent pas être associées une fonction multiplicatrice complexe .
 Utilisez COEF_MULT ou FONC_MULT au lieu de COEF_MULT_C ou FONC_MULT_C.
"""
    ),
}
