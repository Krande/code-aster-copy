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

from ..Utilities import _

cata_msg = {
    1: _(
        """
  Le concept EVOL_CHAR %(k1)s ne contient aucun champ.
"""
    ),
    8: _(
        """
Problème lors du traitement du chargement de type EVOL_CHAR %(k1)s.
L'extraction du chargement de pression a échoué pour l'instant %(r1)f.
Le chargement est mal défini:
- soit %(k1)s n'est pas indexé par l'instant;
- soit le chargement n'a pas été trouvé pour cet instant;
"""
    ),
    12: _(
        """
Problème lors du traitement du chargement de type EVOL_CHAR %(k1)s.
L'extraction du chargement a échoué pour l'instant %(r1)f.
Le chargement est mal défini:
- soit %(k1)s n'est pas indexé par l'instant;
- soit le chargement n'a pas été trouvé pour cet instant;
- soit il manque un champ nécessaire :
    - COEF_H et T_EXT pour ECHANGE
    - FLUN pour FLUX_REP
"""
    ),
}
