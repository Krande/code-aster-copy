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

from ..Objects import (
    PhysicalProblem,
    DiscreteComputation,
)

from ..Utilities import force_list


def calc_matr_elem_ops(self, **args):
    """Execute the command CALC_MATR_ELEM.

    Arguments:
        **args (dict): User's args.

    Returns:
        ElementaryMatrix: elementary matrix
    """

    # Define problem
    model = args["MODELE"]
    mater = args.get("CHAM_MATER")
    cara = args.get("CARA_ELEM")
    phys_pb = PhysicalProblem(model, mater, cara)

    loads = args.get("CHARGE")
    if loads is not None:
        for load in force_list(loads):
            phys_pb.addLoad(load)

    phys_pb.computeListOfLoads()

    disc_comp = DiscreteComputation(phys_pb)

    time = args["INST"]

    myOption = args["OPTION"]

    fourier = args.get("MODE_FOURIER")

    group_ma = args.get("GROUP_MA")
    if group_ma is None:
        group_ma = []
    else:
        group_ma = force_list(group_ma)

    externVar = None
    if phys_pb.getMaterialField() is not None:
        if phys_pb.getMaterialField().hasExternalStateVariable():
            externVar = disc_comp.createExternalStateVariablesField(
                time)

    matr_elem = None
    if myOption == "RIGI_MECA":
        if args["CALC_ELEM_MODELE"] == "OUI":
            matr_elem = disc_comp.elasticStiffnessMatrix(
                time, fourier, group_ma, externVarField=externVar
            )

            matr_rigi_dual = disc_comp.dualStiffnessMatrix()
            matr_elem.addElementaryTerm(
                matr_rigi_dual.getElementaryTerms())
            matr_elem.build()
        else:
            matr_elem = disc_comp.dualStiffnessMatrix()

    elif myOption == "RIGI_GEOM":
        sief_elga = args.get("SIEF_ELGA")
        strx_elga = args.get("STRX_ELGA")
        displ = args.get("DEPL")

        matr_elem = disc_comp.geometricStiffnessMatrix(
            sief_elga, strx_elga, displ, fourier, group_ma)

    elif myOption == "RIGI_ROTA":
        matr_elem = disc_comp.rotationalStiffnessMatrix(group_ma)

    elif myOption == "RIGI_GYRO":
        matr_elem = disc_comp.gyroscopicStiffnessMatrix(group_ma)

    elif myOption == "MECA_GYRO":
        matr_elem = disc_comp.gyroscopicDampingMatrix(group_ma)

    elif myOption == "MASS_MECA":
        matr_elem = disc_comp.massMatrix(False,
                                         group_ma, externVarField=externVar)

    elif myOption == "MASS_MECA_DIAG":
        matr_elem = disc_comp.massMatrix(True,
                                         group_ma, externVarField=externVar)

    elif myOption == "AMOR_MECA":
        massMatrix = args.get("MASS_MECA")
        stiffnessMatrix = args.get("RIGI_MECA")
        matr_elem = disc_comp.dampingMatrix(
            massMatrix, stiffnessMatrix, group_ma, externVarField=externVar)

    elif myOption == "RIGI_MECA_HYST":
        stiffnessMatrix = args["RIGI_MECA"]
        matr_elem = disc_comp.hystereticStiffnessMatrix(
            stiffnessMatrix, group_ma, externVarField=externVar)
        # This part is a remove because complex elem matr does not support real term
        # remove code in hystereticStiffnessMatrix and add this part later
        # matr_rigi_dual = disc_comp.dualStiffnessMatrix()
        # matr_elem.addElementaryTerm(
        #     matr_rigi_dual.getElementaryTerms())
        # matr_elem.build()

    elif myOption == "MASS_THER":
        matr_elem = disc_comp.linearCapacityMatrix(time, group_ma,
                                                   externVarField=externVar)
        matr_elem *= 1.0/args.get("INCR_INST")

    elif myOption == "RIGI_THER":
        matr_elem = disc_comp.linearConductivityMatrix(time, fourier,
                                                       group_ma,
                                                       externVarField=externVar)

        matr_rigi_dual = disc_comp.dualConductivityMatrix()
        matr_elem.addElementaryTerm(matr_rigi_dual.getElementaryTerms())

        matr_rigi_exch = disc_comp.exchangeThermalMatrix(time)
        matr_elem.addElementaryTerm(matr_rigi_exch.getElementaryTerms())

        matr_elem.build()
    elif myOption == "MASS_ACOU":
        matr_elem = disc_comp.compressibilityMatrix(group_ma)

    elif myOption == "AMOR_ACOU":
        matr_elem = disc_comp.impedanceMatrix()

    elif myOption == "RIGI_ACOU":
        matr_elem = disc_comp.linearMobilityMatrix(group_ma)

        matr_rigi_dual = disc_comp.dualMobilityMatrix()
        matr_elem.addElementaryTerm(
            matr_rigi_dual.getElementaryTerms())
        matr_elem.build()

    elif myOption == "RIGI_FLUI_STRU":
        matr_elem = disc_comp.fluidStrucutreStiffnessMatrix(groupOfCells=group_ma,
                                                            externVarField=externVar)

        matr_rigi_dual = disc_comp.dualStiffnessMatrix()
        matr_elem.addElementaryTerm(
            matr_rigi_dual.getElementaryTerms())
        matr_elem.build()

    elif myOption == "MASS_FLUI_STRU":
        matr_elem = disc_comp.fluidStrucutreMassMatrix(groupOfCells=group_ma,
                                                       externVarField=externVar)

    elif myOption == "IMPE_MECA":
        matr_elem = disc_comp.impedanceBoundaryMatrix(group_ma)

    elif myOption == "ONDE_FLUI":
        matr_elem = disc_comp.impedanceWaveMatrix(group_ma)

    else:
        raise RuntimeError("Option %s not implemented" % (myOption))

    return matr_elem
