# coding: utf-8

# Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

from libaster import deleteTemporaryObjects, setFortranLoggingLevel, resetFortranLoggingLevel

from ..Cata.Syntax import _F
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
from ..Utilities import force_list


def affe_materiau_ops(self, **args):
    """Execute the command.

    Arguments:
        **args (dict): User's keywords.
    """

    args = _F(args)

    setFortranLoggingLevel(args["INFO"])

    # create result
    mesh = None
    model = args.get("MODELE")
    if "MAILLAGE" in args:
        mesh = args["MAILLAGE"]
    else:
        mesh = model.getMesh()
    material = MaterialField(mesh)
    if model is not None:
        material.setModel(model)

    fkw = args["AFFE"]
    if isinstance(fkw, dict):
        _addMaterial(material, fkw)
    elif type(fkw) in (list, tuple):
        for curDict in fkw:
            _addMaterial(material, curDict)
    else:
        raise TypeError(
            "Unexpected type: {0!r} {1}".format(fkw, type(fkw)))

    fkw = args.get("AFFE_COMPOR")
    if fkw is not None:
        if isinstance(fkw, dict):
            _addBehaviour(material, fkw)
        elif type(fkw) in (list, tuple):
            for curDict in fkw:
                _addBehaviour(material, curDict)
        else:
            raise TypeError(
                "Unexpected type: {0!r} {1}".format(fkw, type(fkw)))

    externalVarOnMesh = ListOfExternalStateVariables(mesh)
    fkw = args.get("AFFE_VARC")
    if fkw is not None:
        if isinstance(fkw, dict):
            _addExternalStateVariable(material, externalVarOnMesh, fkw, mesh)
        elif type(fkw) in (list, tuple):
            for curDict in fkw:
                _addExternalStateVariable(material, externalVarOnMesh, curDict, mesh)
        else:
            raise TypeError(
                "Unexpected type: {0!r} {1}".format(fkw, type(fkw)))

    material = MaterialFieldBuilder.build(material, externalVarOnMesh)

    resetFortranLoggingLevel()
    deleteTemporaryObjects()

    return material


def _addBehaviour(material, fkw):
    kwTout = fkw.get("TOUT")
    kwGrMa = fkw.get("GROUP_MA")
    compor = fkw["COMPOR"]

    if kwTout is not None:
        material.addBehaviourOnMesh(compor)
    elif kwGrMa is not None:
        kwGrMa = force_list(kwGrMa)
        for grp in kwGrMa:
            material.addBehaviourOnGroupOfCells(compor, grp)
    else:
        raise TypeError("At least {0} or {1} is required"
                        .format("TOUT", "GROUP_MA"))


def _addExternalStateVariable(material, externalVarOnMesh, fkw, mesh):
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
        material.addDependency(chamGd)

    if evol is not None:
        material.addDependency(evol)
        evolParam = EvolutionParameter(evol)
        nomCham = fkw.get("NOM_CHAM")
        if nomCham is not None:
            evolParam.setFieldName(nomCham)
        foncInst = fkw.get("FONC_INST")
        if foncInst is not None:
            evolParam.setTimeFunction(foncInst)

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
        raise RuntimeError("MAILLE is no more supported")
    elif kwGrMa is not None:
        kwGrMa = force_list(kwGrMa)
        for grp in kwGrMa:
            externalVarOnMesh.addExternalStateVariableOnGroupOfCells(
                externalVar, grp)
    else:
        externalVarOnMesh.addExternalStateVariableOnMesh(externalVar)


def _addMaterial(material, fkw):
    kwTout = fkw.get("TOUT")
    kwGrMa = fkw.get("GROUP_MA")
    kwMail = fkw.get("MAILLE")
    mater = fkw["MATER"]

    if type(mater) is not list:
        mater = list(mater)

    if kwTout is not None:
        material.addMaterialsOnMesh(mater)
    elif kwGrMa is not None:
        kwGrMa = force_list(kwGrMa)
        material.addMaterialsOnGroupOfCells(mater, kwGrMa)
    elif kwMail is not None:
        raise RuntimeError("MAILLE is no more supported")
    else:
        raise TypeError("At least {0} or {1} is required"
                        .format("TOUT", "GROUP_MA"))
