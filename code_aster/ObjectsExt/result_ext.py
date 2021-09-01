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
from libaster import Result, AsterError

from ..Utilities import injector, logger
from ..Objects.Serialization import InternalStateBuilder


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
        # list of ranks
        self._st["rank"] = result.getRanks()
        # list of Model objects
        self._st["model"] = []
        # list of MaterialField objects
        self._st["mater"] = []
        # list of ElementaryCharacteristics
        self._st["cara_elem"] = []
        for i in self._st["rank"]:
            try:
                self._st["model"].append(result.getModel(i))
            except AsterError:
                pass
            try:
                self._st["mater"].append(result.getMaterialField(i))
            except AsterError:
                pass
            try:
                self._st["cara_elem"].append(result.getElementaryCharacteristics(i))
            except AsterError:
                pass

        if len(self._st["rank"]) != len(self._st["model"]):
            logger.debug(
                f"Inconsistent definition of models: "
                f"{len(self._st['rank'])} ranks, {len(self._st['model'])} models"
            )
            self._st["model"] = []
        if len(self._st["rank"]) != len(self._st["mater"]):
            logger.debug(
                f"Inconsistent definition of materials fields: "
                f"{len(self._st['rank'])} ranks, {len(self._st['mater'])} materials"
            )
            self._st["mater"] = []
        if len(self._st["cara_elem"]) > 0 and len(self._st["rank"]) != len(self._st["cara_elem"]):
            logger.debug(
                f"Inconsistent definition of elementary characteristics fields: "
                f"{len(self._st['rank'])} ranks, {len(self._st['cara_elem'])} elementary characteristics"
            )
            self._st["cara_elem"] = []
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
                result.setModel(self._st["model"][i], rank)
            if self._st["mater"]:
                result.setMaterialField(self._st["mater"][i], rank)
            if len(self._st["cara_elem"]) > 0 and self._st["cara_elem"]:
                result.setElementaryCharacteristics(self._st["cara_elem"][i], rank)
        if self._st["model"]:
            result.build()


@injector(Result)
class ExtendedResult:

    cata_sdj = "SD.sd_resultat.sd_resultat"
    internalStateBuilder = ResultStateBuilder

    def LIST_CHAMPS(self):
        return aster.GetResu(self.getName(), "CHAMPS")

    def LIST_VARI_ACCES(self):
        return aster.GetResu(self.getName(), "VARI_ACCES")

    def LIST_PARA(self):
        return aster.GetResu(self.getName(), "PARAMETRES")
