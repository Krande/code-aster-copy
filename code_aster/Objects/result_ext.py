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
:py:class:`Result` --- Results container
**************************************************
"""

import aster
from libaster import MaterialField, Model, Result

from ..Utilities import injector
from .Serialization import InternalStateBuilder


class ResultStateBuilder(InternalStateBuilder):
    """Class that returns the internal state of a *Result*."""

    def save(self, result):
        """Return the internal state of a *Result* to be pickled.

        Arguments:
            result (*Result*): The *Result* object to be pickled.

        Returns:
            *InternalStateBuilder*: The internal state itself.
        """
        super().save(result)
        self._st["model"] = []
        self._st["mater"] = []
        self._st["rank"] = result.getRanks()
        for i in self._st["rank"]:
            try:
                self._st["model"].append(result.getModel(i))
            except RuntimeError:
                pass
            try:
                self._st["mater"].append(result.getMaterialField(i))
            except RuntimeError:
                pass
        if len(self._st["rank"]) != len(self._st["model"]):
            self._st["model"] = []
        if len(self._st["rank"]) != len(self._st["mater"]):
            self._st["mater"] = []
        return self

    def restore(self, result):
        """Restore the *DataStructure* content from the previously saved internal
        state.

        Arguments:
            result (*DataStructure*): The *DataStructure* object to be pickled.
        """
        super().restore(result)
        for i, rank in enumerate(self._st["rank"]):
            if self._st["model"]:
                result.addModel(self._st["model"][i], rank)
            if self._st["mater"]:
                result.addMaterialField(self._st["mater"][i], rank)


@injector(Result)
class ExtendedResult(object):

    cata_sdj = "SD.sd_resultat.sd_resultat"
    internalStateBuilder = ResultStateBuilder

    def LIST_CHAMPS (self) :
        return aster.GetResu(self.getName(), "CHAMPS")

    def LIST_VARI_ACCES (self):
        return aster.GetResu(self.getName(), "VARI_ACCES")

    def LIST_PARA (self):
        return aster.GetResu(self.getName(), "PARAMETRES")
