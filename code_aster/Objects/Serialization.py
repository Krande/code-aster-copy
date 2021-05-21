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

from ..Utilities import logger


class InternalStateBuilder:
    """Class that returns the internal state of a *DataStructure* to be pickled.

    Attributes:
        _st (dict): Internal state storage.
    """

    def __init__(self):
        self._st = {"deps": []}

    def __getstate__(self):
        return self._st

    def __setstate__(self, state):
        self._st = state

    def save(self, obj):
        """Return the internal state of a *DataStructure* to be pickled.

        Arguments:
            obj (*DataStructure*): The *DataStructure* object to be pickled.

        Returns:
            *InternalStateBuilder*: The internal state itself.
        """
        assert hasattr(obj, "getType"), f"not a DataStructure: {obj}"
        self._st["deps"] = obj.getDependencies()
        logger.debug(f"saving dependencies of {obj}: {self._st['deps']}")
        return self

    def restore(self, obj):
        """Restore the *DataStructure* content from the previously saved internal
        state.

        Arguments:
            obj (*DataStructure*): The *DataStructure* object to be pickled.
        """
        logger.debug(f"restoring dependencies for {obj}: {self._st['deps']}")
        for ref in self._st["deps"]:
            obj.addDependency(ref)
