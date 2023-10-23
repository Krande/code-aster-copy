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

from .base_opers_manager import BaseOperatorsManager
from ...Utilities import no_new_attributes


class ThermalOperatorsManager(BaseOperatorsManager):
    """Solve an iteration."""

    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        super().__init__()

    def getResidual(self, scaling=1.0):
        """Computes the functional."""
        raise NotImplementedError

    def getJacobian(self, matrix_type):
        """Computes the jacobian."""
        raise NotImplementedError
