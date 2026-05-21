# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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

"""
Fonctionnalités utilitaires pour CALC_ENDO
"""

from ...Utilities import logger
from ...Messages import UTMESS

from ...CodeCommands import CALC_TABLE


def get_obs_values(_ctrl_resu, _name_obs):
    """Get values of a given OBSERVATION

    Args:
        _ctrl_resu (*Table*): observation table
        _name_obs (str): name of the quantity to extract from the table

    Returns:
        _ctrl_obs_values (list): list of obervation table values
    """

    _ctrl_obs = CALC_TABLE(
        TABLE=_ctrl_resu,
        ACTION=(
            _F(OPERATION="FILTRE", NOM_PARA="NOM_OBSERVATION", VALE_K=_name_obs),
            _F(OPERATION="EXTR", NOM_PARA=("INST", "VALE")),
        ),
    )

    _ctrl_obs_values = _ctrl_obs.EXTR_TABLE().values()["VALE"]

    return _ctrl_obs_values
