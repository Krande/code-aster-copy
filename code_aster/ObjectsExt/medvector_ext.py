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

# person_in_charge: mathieu.courtois@edf.fr
"""
:py:class:`DataStructure` --- Base of all objects
*************************************************
"""

import libaster
from libaster import MedVector

from ..Utilities import injector
from ..Objects.Serialization import InternalStateBuilder


@injector(MedVector)
class ExtendedMedVector:
    """This class defines the base class of the MedVector."""

    def __getstate__(self):
        """Return internal state.

        Derivated *DataStructure* types should defined a dedicated *InternalStateBuilder*
        class to serialize its specific content.
        """
        return (
            self.getValues(),
            self.getCumulatedSizesVector(),
            self.getComponentVector(),
            self.getComponentNumber(),
            self.getComponentName(),
        )

    def __setstate__(self, state):
        """Restore internal state.

        Arguments:
            state (*InternalStateBuilder*): Internal state.
        """
        self.__init__()
        self.setValues(state[0])
        self.setCumulatedSizesVector(state[1])
        self.setComponentVector(state[2])
        self.setComponentNumber(state[3])
        self.setComponentName(state[4])
