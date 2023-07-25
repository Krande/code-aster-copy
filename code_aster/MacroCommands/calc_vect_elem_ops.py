# coding: utf-8

# Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

from ..Objects import (
    PhysicalProblem,
    DiscreteComputation,
    ElementaryVectorDisplacementReal,
    ElementaryVectorTemperatureReal,
    ElementaryVectorPressureComplex,
    FieldOnNodesReal,
)

from ..Utilities import force_list


def calc_vect_elem_ops(self, **args):
    """Execute the command CALC_VECT_ELEM.

    Arguments:
        **args (dict): User's args.

    Returns:
        ElementaryVector: elementary vector
    """

    loads = force_list(args.get("CHARGE"))
    myOption = args["OPTION"]

    fourier = args.get("MODE_FOURIER")
    if not fourier:
        fourier = 0
    time = args.get("INST")

    model = args.get("MODELE")
    if model is None:
        model = loads[0].getModel()

    assert model is not None

    mater = args.get("CHAM_MATER")
    cara = args.get("CARA_ELEM")

    # Define problem
    phys_pb = PhysicalProblem(model, mater, cara)

    for load in loads:
        phys_pb.addLoad(load)

    phys_pb.computeListOfLoads()

    disc_comp = DiscreteComputation(phys_pb)

    varc = None
    if mater and mater.hasExternalStateVariable():
        varc = phys_pb.getExternalStateVariables(time)

    if myOption == "CHAR_MECA":
        vect_elem = ElementaryVectorDisplacementReal(
            phys_pb.getModel(),
            phys_pb.getMaterialField(),
            phys_pb.getElementaryCharacteristics(),
            phys_pb.getListOfLoads(),
        )
    elif myOption == "CHAR_THER":
        vect_elem = ElementaryVectorTemperatureReal(
            phys_pb.getModel(),
            phys_pb.getMaterialField(),
            phys_pb.getElementaryCharacteristics(),
            phys_pb.getListOfLoads(),
        )
    elif myOption == "CHAR_ACOU":
        vect_elem = ElementaryVectorPressureComplex(
            phys_pb.getModel(),
            phys_pb.getMaterialField(),
            phys_pb.getElementaryCharacteristics(),
            phys_pb.getListOfLoads(),
        )
    else:
        raise RuntimeError("Option %s not implemented" % (myOption))

    vect_elem.prepareCompute(myOption)

    neum_elem = disc_comp.getNeumannForces(time, mode=fourier, varc_curr=varc, assembly=False)
    vect_elem.addElementaryTerm(neum_elem.getElementaryTerms())

    dual_elem = disc_comp.getImposedDualBC(time, assembly=False)
    vect_elem.addElementaryTerm(dual_elem.getElementaryTerms())

    if phys_pb.getModel().isThermal():
        exch_elem = disc_comp.getThermalExchangeForces(
            FieldOnNodesReal(phys_pb.getModel()), time, assembly=False
        )
        vect_elem.addElementaryTerm(exch_elem.getElementaryTerms())

    if "SOUS_STRUC" in args:
        load_struc = {}
        for struc in force_list(args["SOUS_STRUC"]):
            super_cells = []
            if "SUPER_MAILLE" in struc:
                super_cells = force_list(struc["SUPER_MAILLE"])
            load_struc[struc["CAS_CHARGE"]] = super_cells
        vect_elem.addSubstructuring(load_struc)

    vect_elem.build()

    return vect_elem
