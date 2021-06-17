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

#
# This file is part of Code_Aster.
#
# Code_Aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Code_Aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.

# person_in_charge: nicolas.sellenet@edf.fr

from ..Objects import (ConcreteDryingExternalStateVariable,
                       ConcreteHydratationExternalStateVariable,
                       CorrosionExternalStateVariable, EvolutionParameter,
                       GeometryExternalStateVariable,
                       ListOfExternalStateVariables, IrradiationExternalStateVariable,
                       IrreversibleDeformationExternalStateVariable, MaterialField,
                       MaterialFieldBuilder, Neutral1ExternalStateVariable,
                       Neutral2ExternalStateVariable, Neutral3ExternalStateVariable,
                       SteelPhasesExternalStateVariable, TemperatureExternalStateVariable,
                       TotalFluidPressureExternalStateVariable,
                       VolumetricDeformationExternalStateVariable,
                       ZircaloyPhasesExternalStateVariable)
from ..Supervis import ExecuteCommand
from ..Utilities import force_list


class MaterialAssignment(ExecuteCommand):
    """Assign the :class:`~code_aster.Objects.Material` properties on the
    :class:`~code_aster.Objects.Mesh` that creates a
    :class:`~code_aster.Objects.MaterialField` object.
    """
    command_name = "AFFE_MATERIAU"

    def create_result(self, keywords):
        """Initialize the :class:`~code_aster.Objects.Mesh`.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        mesh = None
        model = keywords.get("MODELE")
        if "MAILLAGE" in keywords:
            mesh = keywords["MAILLAGE"]
        else:
            mesh = model.getMesh()
        self._result = MaterialField(mesh)
        if model is not None:
            self._result.setModel(model)

    def exec_(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        fkw = keywords["AFFE"]
        if isinstance(fkw, dict):
            self._addMaterial(fkw)
        elif type(fkw) in (list, tuple):
            for curDict in fkw:
                self._addMaterial(curDict)
        else:
            raise TypeError("Unexpected type: {0!r} {1}".format(fkw, type(fkw)))

        fkw = keywords.get("AFFE_COMPOR")
        if fkw is not None:
            if isinstance(fkw, dict):
                self._addBehaviour(fkw)
            elif type(fkw) in (list, tuple):
                for curDict in fkw:
                    self._addBehaviour(curDict)
            else:
                raise TypeError("Unexpected type: {0!r} {1}".format(fkw, type(fkw)))

        mesh = None
        if "MAILLAGE" in keywords:
            mesh = keywords["MAILLAGE"]
        else:
            mesh = keywords["MODELE"].getMesh()

        externalVarOnMesh = ListOfExternalStateVariables(mesh)
        fkw = keywords.get("AFFE_VARC")
        if fkw is not None:
            if isinstance(fkw, dict):
                self._addExternalStateVariable(externalVarOnMesh, fkw, mesh)
            elif type(fkw) in (list, tuple):
                for curDict in fkw:
                    self._addExternalStateVariable(externalVarOnMesh, curDict, mesh)
            else:
                raise TypeError("Unexpected type: {0!r} {1}".format(fkw, type(fkw)))

        self._result = MaterialFieldBuilder.build(self._result, externalVarOnMesh)

    def _addBehaviour(self, fkw):
        kwTout = fkw.get("TOUT")
        kwGrMa = fkw.get("GROUP_MA")
        mater = fkw["COMPOR"]

        if kwTout is not None:
            self._result.addBehaviourOnMesh(mater)
        elif kwGrMa is not None:
            kwGrMa = force_list(kwGrMa)
            for grp in kwGrMa:
                self._result.addBehaviourOnGroupOfCells(mater, grp)
        else:
            raise TypeError("At least {0} or {1} is required"
                            .format("TOUT", "GROUP_MA"))

    def _addExternalStateVariable(self, externalVarOnMesh, fkw, mesh):
        kwTout = fkw.get("TOUT")
        kwGrMa = fkw.get("GROUP_MA")
        kwMail = fkw.get("MAILLE")
        nomVarc = fkw["NOM_VARC"]
        chamGd = fkw.get("CHAM_GD")
        valeRef = fkw.get("VALE_REF")
        evol = fkw.get("EVOL")

        obj = None
        if nomVarc == "TEMP":
            obj = TemperatureExternalStateVariable
        elif nomVarc == "GEOM":
            obj = GeometryExternalStateVariable
        elif nomVarc == "CORR":
            obj = CorrosionExternalStateVariable
        elif nomVarc == "EPSA":
            obj = IrreversibleDeformationExternalStateVariable
        elif nomVarc == "HYDR":
            obj = ConcreteHydratationExternalStateVariable
        elif nomVarc == "IRRA":
            obj = IrradiationExternalStateVariable
        elif nomVarc == "M_ACIER":
            obj = SteelPhasesExternalStateVariable
        elif nomVarc == "M_ZIRC":
            obj = ZircaloyPhasesExternalStateVariable
        elif nomVarc == "NEUT1":
            obj = Neutral1ExternalStateVariable
        elif nomVarc == "NEUT2":
            obj = Neutral2ExternalStateVariable
        elif nomVarc == "NEUT3":
            obj = Neutral3ExternalStateVariable
        elif nomVarc == "SECH":
            obj = ConcreteDryingExternalStateVariable
        elif nomVarc == "PTOT":
            obj = TotalFluidPressureExternalStateVariable
        elif nomVarc == "DIVU":
            obj = VolumetricDeformationExternalStateVariable
        else:
            raise TypeError("Input Variable not allowed")

        externalVar = obj(mesh)
        if valeRef is not None:
            externalVar.setReferenceValue(valeRef)

        if chamGd is not None:
            externalVar.setValue(chamGd)

        if evol is not None:
            evolParam = EvolutionParameter(evol)
            nomCham = fkw.get("NOM_CHAM")
            if nomCham is not None: evolParam.setFieldName(nomCham)
            foncInst = fkw.get("FONC_INST")
            if foncInst is not None: evolParam.setTimeFunction(foncInst)

            prolDroite = fkw.get("PROL_DROITE")
            if prolDroite is not None:
                if prolDroite == "EXCLU":
                    evolParam.prohibitRightExtension()
                if prolDroite == "CONSTANT":
                    evolParam.setConstantRightExtension()
                if prolDroite == "LINEAIRE":
                    evolParam.setLinearRightExtension()

            prolGauche = fkw.get("PROL_GAUCHE")
            if prolGauche is not None:
                if prolGauche == "EXCLU":
                    evolParam.prohibitLeftExtension()
                if prolGauche == "CONSTANT":
                    evolParam.setConstantLeftExtension()
                if prolGauche == "LINEAIRE":
                    evolParam.setLinearLeftExtension()

            externalVar.setEvolutionParameter(evolParam)

        if kwTout is not None:
            externalVarOnMesh.addExternalStateVariableOnMesh(externalVar)
        elif kwMail is not None:
            kwMail = force_list(kwMail)
            for elem in kwMail:
                externalVarOnMesh.addExternalStateVariableOnCell(externalVar, elem)
        elif kwGrMa is not None:
            kwGrMa = force_list(kwGrMa)
            for grp in kwGrMa:
                externalVarOnMesh.addExternalStateVariableOnGroupOfCells(externalVar, grp)
        else:
            externalVarOnMesh.addExternalStateVariableOnMesh(externalVar)

    def _addMaterial(self, fkw):
        kwTout = fkw.get("TOUT")
        kwGrMa = fkw.get("GROUP_MA")
        kwMail = fkw.get("MAILLE")
        mater = fkw["MATER"]
        if type(mater) is not list:
            mater = list(mater)

        if kwTout is not None:
            self._result.addMaterialsOnMesh(mater)
        elif kwGrMa is not None:
            kwGrMa = force_list(kwGrMa)
            self._result.addMaterialsOnGroupOfCells(mater, kwGrMa)
        elif kwMail is not None:
            kwMail = force_list(kwMail)
            self._result.addMaterialsOnCell(mater, kwMail)
        else:
            raise TypeError("At least {0}, {1} or {2} is required"
                            .format("TOUT", "GROUP_MA", "MAILLE"))


AFFE_MATERIAU = MaterialAssignment.run
