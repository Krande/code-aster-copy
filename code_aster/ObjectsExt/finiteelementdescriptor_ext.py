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

"""
:py:class:`FiniteElementDescriptor` --- Finite Element definition
*****************************************
"""

from libaster import FiniteElementDescriptor

from ..Utilities import injector
from ..Objects.Serialization import InternalStateBuilder


class FEDStateBuilder(InternalStateBuilder):
    """Class that returns the internal state of a *FiniteElementDescriptor*."""

    def save(self, fed):
        """Return the internal state of a *FiniteElementDescriptor* to be pickled.

        Arguments:
            fed (*FiniteElementDescriptor*): The *FiniteElementDescriptor* object to be pickled.

        Returns:
            *InternalStateBuilder*: The internal state itself.
        """
        super().save(fed)
        self._st["model"] = fed.getModel()
        return self

    def restore(self, fed):
        """Restore the *DataStructure* content from the previously saved internal
        state.

        Arguments:
            field (*DataStructure*): The *DataStructure* object to be restored.
        """
        super().restore(fed)
        if self._st["model"]:
            fed.setModel(self._st["model"])


@injector(FiniteElementDescriptor)
class ExtendedFiniteElementDescriptor:
    cata_sdj = "SD.sd_ligrel.sd_ligrel"
    internalStateBuilder = FEDStateBuilder

    def __getinitargs__(self):
        """Returns the argument required to reinitialize a
        FiniteElementDescriptor object during unpickling.
        """
        return (self.getName(), self.getMesh())
