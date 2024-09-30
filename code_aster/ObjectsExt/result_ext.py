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
from ..Messages import UTMESS
from ..Utilities import injector, logger, SearchList, is_number, medcoupling as medc
from ..Objects.Serialization import InternalStateBuilder


class ResultStateBuilder(InternalStateBuilder):
    """Class that returns the internal state of a *Result*."""

    def _addFields(self, result, fieldNames, indexes):
        for i in indexes:
            for fieldName in fieldNames:
                try:
                    curField = result.getField(fieldName, i)
                    self._st["fields"][i][fieldName] = curField
                except:
                    pass

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

        indexes = self._st["index"]
        self._st["fields"] = {}
        for i in indexes:
            self._st["fields"][i] = {}
        self._addFields(result, result.getFieldsOnNodesRealNames(), indexes)
        self._addFields(result, result.getFieldsOnNodesComplexNames(), indexes)
        self._addFields(result, result.getFieldsOnCellsRealNames(), indexes)
        self._addFields(result, result.getFieldsOnCellsComplexNames(), indexes)
        self._addFields(result, result.getFieldsOnCellsLongNames(), indexes)
        self._addFields(result, result.getConstantFieldsOnCellsRealNames(), indexes)
        self._addFields(result, result.getConstantFieldsOnCellsChar16Names(), indexes)

        return self

    def restore(self, result):
        """Restore the *DataStructure* content from the previously saved internal
        state.

        Arguments:
            result (*DataStructure*): The *DataStructure* object to be restored.
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

        for fed in self._st["feds"]:
            result.addFiniteElementDescriptor(fed)
        for fnd in self._st["fnds"]:
            result.addEquationNumbering(fnd)
        for index in self._st["fields"]:
            fields = self._st["fields"][index]
            for fieldName in fields:
                result.setField(fields[fieldName], fieldName, index)


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

    def _createIndexFromParameter(self, para, value, crit, prec):
        """
        Create the index corresponding to a given value of an access parameter.

        Arguments:
            para (str) : name of the access parameter (NUME_ORDRE, INST, etc..)
            value (float|int|str) : value of the access parameter
            crit (str) : search criterion ABSOLU or RELATIF
            prec (float) : precision for the search criterion

        Returns:
            index (int) : the corresponding index (index)

        """
        acpara = self.getAccessParameters()
        if para not in acpara:
            UTMESS("F", "RESULT1_8")
        if isinstance(value, int):
            storageIndex = self.createIndexFromParameter(para, value)
        elif isinstance(value, str):
            storageIndex = self.createIndexFromParameter(para, value)
        else:
            raise ValueError(f"Type of access to result is invalid {value!r}")

        return storageIndex

    def _getIndexFromParameter(self, para, value, crit, prec, throw_except):
        """
        Get the index corresponding to a given value of an access parameter.

        Arguments:
            para (str) : name of the access parameter (NUME_ORDRE, INST, etc..)
            value (float|int|str) : value of the access parameter
            crit (str) : search criterion ABSOLU or RELATIF
            prec (float) : precision for the search criterion

        Returns:
            index (int) : the corresponding index (index)

        """

        acpara = self.getAccessParameters()
        if not para in acpara:
            msg = "Missing parameter {}".format(para)
            raise ValueError(msg)

        if para not in acpara:
            UTMESS("F", "RESULT1_8")
        if is_number(value):
            slist = SearchList(acpara[para], prec, crit)
            internalStorage = slist.index(value)

        elif isinstance(value, str):
            slist = acpara[para]
            if throw_except:
                internalStorage = slist.index(value)
            else:
                try:
                    internalStorage = slist.index(value)
                except ValueError:
                    internalStorage = -1
                    return internalStorage
        else:
            raise ValueError(f"Type of access to result is invalid {value!r}")

        return acpara["NUME_ORDRE"][internalStorage]

    def getField(self, name, value=None, para="NUME_ORDRE", crit="RELATIF", prec=1.0e-6):
        """Get the specified field. This is an overlay to existing methods
        for each type of field.

        Arguments:
            name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
            value (float|int|str) : value of the access parameter
            para (str) : name of the access parameter (NUME_ORDRE, INST, etc..)
            crit (str) : search criterion ABSOLU or RELATIF
            prec (float) : precision for the search criterion

        Returns:
            Field***: field to get whit type in (FieldOnNodes***/FieldOnCells***/
            ConstantFieldOnCell***)
        """

        assert crit in ("ABSOLU", "RELATIF")

        if para in ("NUME_ORDRE",):
            storageIndex = value
        else:
            storageIndex = self._getIndexFromParameter(para, value, crit, prec, throw_except=True)

        if storageIndex == -1:
            UTMESS("F", "RESULT1_9")

        names = self.getFieldsOnNodesRealNames()
        if name in names:
            return self._getFieldOnNodesReal(name, storageIndex)

        names = self.getFieldsOnNodesComplexNames()
        if name in names:
            return self._getFieldOnNodesComplex(name, storageIndex)

        names = self.getFieldsOnCellsRealNames()
        if name in names:
            return self._getFieldOnCellsReal(name, storageIndex)

        names = self.getFieldsOnCellsComplexNames()
        if name in names:
            return self._getFieldOnCellsComplex(name, storageIndex)

        names = self.getFieldsOnCellsLongNames()
        if name in names:
            return self._getFieldOnCellsLong(name, storageIndex)

        names = self.getConstantFieldsOnCellsRealNames()
        if name in names:
            return self._getConstantFieldOnCellsReal(name, storageIndex)

        names = self.getConstantFieldsOnCellsChar16Names()
        if name in names:
            return self._getConstantFieldOnCellsChar16(name, storageIndex)

        names = self.getGeneralizedVectorRealNames()
        if name in names:
            return self.getGeneralizedVectorReal(name, storageIndex)

        names = self.getGeneralizedVectorComplexNames()
        if name in names:
            return self.getGeneralizedVectorComplex(name, storageIndex)

        raise KeyError("name of field %s not found" % name)

    def setField(self, field, name, value=None, para="NUME_ORDRE", crit="RELATIF", prec=1.0e-6):
        """Set the specified field. This is an overlay to existing methods
        for each type of field.

        Arguments:
            name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
            field : field
            value (float|int|str) : value of the access parameter
            para (str) : name of the access parameter (NUME_ORDRE, INST, etc..)
            crit (str) : search criterion ABSOLU or RELATIF
            prec (float) : precision for the search criterion

        Returns:
            Nothing
        """

        assert crit in ("ABSOLU", "RELATIF")

        if para in ("NUME_ORDRE",):
            storageIndex = value
        else:
            storageIndex = self._getIndexFromParameter(para, value, crit, prec, throw_except=False)

        if storageIndex < 0:
            storageIndex = self._createIndexFromParameter(para, value, crit, prec)
            if storageIndex < 0:
                raise KeyError("Echec lors de la création du paramètre")

        self._setField(field, name, storageIndex)

    def createMedCouplingResult(self, medmesh=None):
        """Export the result to a new MED container.

        The export is limited to fields on nodes (Real) only.

        Arguments:
            medmesh, optional (*MEDFileUMesh*): The medcoupling support mesh.

        Returns:
            field ( MEDFileData ) : The result in med format ( medcoupling ).
        """

        # Get NUME_ORDRE
        para = self.getAccessParameters()

        ranks = para.get("NUME_ORDRE")
        assert ranks is not None

        # Get the variable to be associated to time in the med file
        times = None
        for i in ("INST", "FREQ", "NUME_MODE"):
            if i in para:
                times = para.get(i)
                break
        assert times is not None

        # Only works for fields on nodes real
        save_fields = self.getFieldsOnNodesRealNames()
        if len(save_fields) == 0:
            msg = "None of the fields can be exported to medcoupling"
            raise RuntimeError(msg)

        # Init medcoupling objets
        medresult = medc.MEDFileData()
        meshes = medc.MEDFileMeshes()
        fields = medc.MEDFileFields()

        # Set mesh
        medmesh = medmesh or self.getMesh().createMedCouplingMesh()
        meshes.pushMesh(medmesh)

        for fname in save_fields:
            fmts = medc.MEDFileFieldMultiTS()
            for rank, time in zip(ranks, times):
                medcfield = self.getField(fname, rank).createMedCouplingField(medmesh)
                medcfield.setTime(rank, 0, time)
                fmts.pushBackTimeStep(medcfield)
            fields.pushField(fmts)

        medresult.setMeshes(meshes)
        medresult.setFields(fields)
        return medresult
