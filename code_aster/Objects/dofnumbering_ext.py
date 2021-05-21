# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
:py:class:`DOFNumbering` --- DOFNumbering definition
*****************************************
"""

import aster
from libaster import DOFNumbering

from ..Utilities import injector
from .Serialization import InternalStateBuilder


class DOFNumberingStateBuilder(InternalStateBuilder):
    """Class that returns the internal state of a *DOFNumbering* to be pickled."""

    def save(self, numbering):
        """Return the internal state of a *DOFNumbering* to be pickled.

        Arguments:
            numbering (*DOFNumbering*): The *DOFNumbering* object to be pickled.

        Returns:
            *InternalStateBuilder*: The internal state itself.
        """
        super().save(numbering)
        self._st["model"] = numbering.getModel()
        return self

    def restore(self, numbering):
        """Restore the *DOFNumbering* content from the previously saved internal
        state.

        Arguments:
            numbering (*DOFNumbering*): The *DataStructure* object to be pickled.
        """
        super().restore(numbering)
        numbering.setModel(self._st["model"])


@injector(DOFNumbering)
class ExtendedDOFNumbering(object):
    cata_sdj = "SD.sd_nume_ddl.sd_nume_ddl"
    internalStateBuilder = DOFNumberingStateBuilder
