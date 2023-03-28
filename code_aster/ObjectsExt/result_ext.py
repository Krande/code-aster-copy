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
:py:class:`Result` --- Results container
**************************************************
"""

import aster
from libaster import Result

from ..Utilities import injector, logger, SearchList
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
        # mesh
        self._st["mesh"] = result.getMesh()
        # list of indexs
        self._st["index"] = result.getIndexes()
        # list of Model objects
        self._st["model"] = []
        # list of MaterialField objects
        self._st["mater"] = []
        # list of ElementaryCharacteristics
        self._st["cara_elem"] = []
        # list of list of loads
        self._st["loads"] = []
        # list of ligrel
        self._st["feds"] = result.getFiniteElementDescriptors()
        # list of nume_equa
        self._st["fnds"] = result.getEquationNumberings()
        for i in self._st["index"]:
            if result.hasModel(i):
                self._st["model"].append(result.getModel(i))
            if result.hasMaterialField(i):
                self._st["mater"].append(result.getMaterialField(i))
            if result.hasElementaryCharacteristics(i):
                self._st["cara_elem"].append(result.getElementaryCharacteristics(i))
            if result.hasListOfLoads(i):
                self._st["loads"].append(result.getListOfLoads(i))

        if len(self._st["index"]) != len(self._st["model"]):
            logger.debug(
                f"Inconsistent definition of models: "
                f"{len(self._st['index'])} indexs, {len(self._st['model'])} models"
            )
            self._st["model"] = []
        if len(self._st["index"]) != len(self._st["mater"]):
            logger.debug(
                f"Inconsistent definition of materials fields: "
                f"{len(self._st['index'])} indexs, {len(self._st['mater'])} materials"
            )
            self._st["mater"] = []
        if len(self._st["cara_elem"]) > 0 and len(self._st["index"]) != len(self._st["cara_elem"]):
            logger.debug(
                f"Inconsistent definition of elementary characteristics fields: "
                f"{len(self._st['index'])} indexs, {len(self._st['cara_elem'])} elementary characteristics"
            )
            self._st["cara_elem"] = []
        if len(self._st["loads"]) > 0 and len(self._st["index"]) != len(self._st["loads"]):
            logger.debug(
                f"Inconsistent definition of list of loads: "
                f"{len(self._st['index'])} indexs, {len(self._st['loads'])} list of loads"
            )
            self._st["loads"] = []
        return self

    def restore(self, result):
        """Restore the *DataStructure* content from the previously saved internal
        state.

        Arguments:
            result (*DataStructure*): The *DataStructure* object to be pickled.
        """
        super().restore(result)
        result.setMesh(self._st["mesh"])
        for i, index in enumerate(self._st["index"]):
            if self._st["model"]:
                result.setModel(self._st["model"][i], index)
            if self._st["mater"]:
                result.setMaterialField(self._st["mater"][i], index)
            if len(self._st["cara_elem"]) > 0 and self._st["cara_elem"]:
                result.setElementaryCharacteristics(self._st["cara_elem"][i], index)
            if len(self._st["loads"]) > 0 and self._st["loads"]:
                result.setListOfLoads(self._st["loads"][i], index)

        if result.getMesh():
            result.build(self._st["feds"], self._st["fnds"])


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

    def _getIndexFromParameter(self, para, value, crit, prec):
        """
        Get the index corresponding to a given value of an access parameter.

        Arguments:
            para (str) : name of the access parameter (NUME_ORDRE, INST, etc..)
            value (float|int) : value of the access parameter
            crit (str) : search criterion ABSOLU or RELATIF
            prec (float) : precision for the search criterion

        Returns:
            index (int) : the corresponding index (index)

        """
        acpara = self.getAccessParameters()

        if not para in acpara:
            msg = "Missing parameter {}".format(para)
            raise ValueError(msg)

        slist = SearchList(acpara[para], prec, crit)
        idx = slist.index(value)

        return acpara["NUME_ORDRE"][idx]

    def getField(self, name, value=None, para="NUME_ORDRE", crit="RELATIF", prec=1.0e-6):
        """Get the specified field. This is an overlay to existing methods
        for each type of field.

        Arguments:
            name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
            value (float|int) : value of the access parameter
            para (str) : name of the access parameter (NUME_ORDRE, INST, etc..)
            crit (str) : search criterion ABSOLU or RELATIF
            prec (float) : precision for the search criterion

        Returns:
            Field***: field to get whit type in (FieldOnNodes***/FieldOnCells***/
            ConstantFieldOnCell***)
        """

        assert crit in ("ABSOLU", "RELATIF")

        if para in ("NUME_ORDRE"):
            index = value
        else:
            index = self._getIndexFromParameter(para, value, crit, prec)

        names = self.getFieldsOnNodesRealNames()
        if name in names:
            return self.getFieldOnNodesReal(name, index)

        names = self.getFieldsOnNodesComplexNames()
        if name in names:
            return self.getFieldOnNodesComplex(name, index)

        names = self.getFieldsOnCellsRealNames()
        if name in names:
            return self.getFieldOnCellsReal(name, index)

        names = self.getFieldsOnCellsComplexNames()
        if name in names:
            return self.getFieldOnCellsComplex(name, index)

        names = self.getFieldsOnCellsLongNames()
        if name in names:
            return self.getFieldOnCellsLong(name, index)

        names = self.getConstantFieldsOnCellsRealNames()
        if name in names:
            return self.getConstantFieldOnCellsReal(name, index)

        names = self.getConstantFieldsOnCellsChar16Names()
        if name in names:
            return self.getConstantFieldOnCellsChar16(name, index)

        names = self.getGeneralizedVectorRealNames()
        if name in names:
            return self.getGeneralizedVectorReal(name, index)

        names = self.getGeneralizedVectorComplexNames()
        if name in names:
            return self.getGeneralizedVectorComplex(name, index)

        raise KeyError("name of field %s not found" % name)
