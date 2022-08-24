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

# person_in_charge: nicolas.sellenet@edf.fr

from ..Objects import (
    ElementaryMatrixDisplacementComplex,
    ElementaryMatrixDisplacementReal,
    ElementaryMatrixPressureComplex,
    ElementaryMatrixTemperatureReal,
    PhysicalProblem,
    MechanicalLoadComplex,
    DiscreteComputation,
)
from ..Supervis import ExecuteCommand
from ..Utilities import force_list


class ComputeElementaryMatrix(ExecuteCommand):

    """Command that creates evolutive results."""

    command_name = "CALC_MATR_ELEM"

    def use_cpp(self, keywords):
        """Use or not new c++ commands"""

        # Case not yet supported - to fix

        myOption = keywords["OPTION"]
        if myOption not in ("RIGI_MECA", "MASS_MECA", "AMOR_MECA",
                            "RIGI_THER", "MASS_THER", "RIGI_MECA_HYST",
                            "MASS_ACOU", "RIGI_ACOU", "AMOR_ACOU",):
            return False

        strucElem = keywords.get("CARA_ELEM")
        # AMOR_MECA in C++ doesnot work for some structural elements
        if strucElem is not None:
            if myOption in ("AMOR_MECA"):
                return False

        loads = keywords.get("CHARGE")
        if loads is not None:
            for load in force_list(loads):
                if isinstance(load, MechanicalLoadComplex):
                    return False

        return True

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        myOption = keywords["OPTION"]
        if myOption in (
            "AMOR_MECA",
            "AMOR_MECA_ABSO",
            "IMPE_MECA",
            "MASS_FLUI_STRU",
            "MASS_MECA",
            "MASS_MECA_DIAG",
            "MECA_GYRO",
            "ONDE_FLUI",
            "RIGI_FLUI_STRU",
            "RIGI_GEOM",
            "RIGI_GYRO",
            "RIGI_MECA",
            "RIGI_ROTA",
        ):
            self._result = ElementaryMatrixDisplacementReal()
        elif myOption == "RIGI_MECA_HYST":
            self._result = ElementaryMatrixDisplacementComplex()
        elif myOption in ("RIGI_THER", "MASS_THER"):
            self._result = ElementaryMatrixTemperatureReal()
        elif myOption in ("RIGI_ACOU", "MASS_ACOU", "AMOR_ACOU"):
            self._result = ElementaryMatrixPressureComplex()

    def exec_(self, keywords):
        """Override default _exec in case of some options"""
        # Compute option
        if self.use_cpp(keywords):
            # Define problem
            model = keywords["MODELE"]
            mater = keywords.get("CHAM_MATER")
            cara = keywords.get("CARA_ELEM")
            phys_pb = PhysicalProblem(model, mater, cara)

            loads = keywords.get("CHARGE")
            if loads is not None:
                for load in force_list(loads):
                    phys_pb.addLoad(load)

            phys_pb.computeListOfLoads()

            disc_comp = DiscreteComputation(phys_pb)

            time = keywords["INST"]
            delta_time = keywords.get("INCR_INST")
            myOption = keywords["OPTION"]

            fourier = keywords.get("MODE_FOURIER")

            group_ma = keywords.get("GROUP_MA")
            if group_ma is None:
                group_ma = []
            else:
                group_ma = force_list(group_ma)

            externVar = None
            if phys_pb.getMaterialField() is not None:
                if phys_pb.getMaterialField().hasExternalStateVariable():
                    externVar = disc_comp.createExternalStateVariablesField(
                        time)

            if myOption == "RIGI_MECA":
                if keywords["CALC_ELEM_MODELE"] == "OUI":
                    self._result = disc_comp.elasticStiffnessMatrix(
                        time, fourier, group_ma, externVarField=externVar
                    )

                    matr_rigi_dual = disc_comp.dualStiffnessMatrix()
                    self._result.addElementaryTerm(
                        matr_rigi_dual.getElementaryTerms())
                    self._result.build()
                else:
                    self._result = disc_comp.dualStiffnessMatrix()

            elif myOption == "MASS_MECA":
                self._result = disc_comp.massMatrix(
                    group_ma, externVarField=externVar)

            elif myOption == "AMOR_MECA":
                massMatrix = keywords.get("MASS_MECA")
                stiffnessMatrix = keywords.get("RIGI_MECA")
                self._result = disc_comp.dampingMatrix(
                    massMatrix, stiffnessMatrix, group_ma, externVarField=externVar)

            elif myOption == "RIGI_MECA_HYST":
                stiffnessMatrix = keywords["RIGI_MECA"]
                self._result = disc_comp.complexStiffnessMatrix(
                    stiffnessMatrix, group_ma, externVarField=externVar)

            elif myOption == "MASS_THER":
                self._result = disc_comp.linearCapacityMatrix(time, delta_time, 1.0,
                                                              group_ma,
                                                              externVarField=externVar)

            elif myOption == "RIGI_THER":
                self._result = disc_comp.linearConductivityMatrix(time, delta_time, 1.0,
                                                                  fourier,
                                                                  group_ma,
                                                                  externVarField=externVar)
            elif myOption == "MASS_ACOU":
                self._result = disc_comp.compressibilityMatrix(group_ma)

            elif myOption == "AMOR_ACOU":
                self._result = disc_comp.impedanceMatrix()

            elif myOption == "RIGI_ACOU":
                self._result = disc_comp.linearMobilityMatrix(group_ma)

                matr_rigi_dual = disc_comp.dualMobilityMatrix()
                self._result.addElementaryTerm(
                    matr_rigi_dual.getElementaryTerms())
                self._result.build()

            else:
                raise RuntimeError("Option %s not implemented" % (myOption))
        else:
            super(ComputeElementaryMatrix, self).exec_(keywords)

    def post_exec(self, keywords):
        """Store references to ElementaryMatrix objects.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords, changed
                in place.
        """

        if not self.use_cpp(keywords):
            self._result.setModel(keywords["MODELE"])

            chamMater = keywords.get("CHAM_MATER")
            if chamMater is not None:
                self._result.setMaterialField(chamMater)

            caraElem = keywords.get("CARA_ELEM")
            if caraElem is not None:
                self._result.setElementaryCharacteristics(caraElem)
            self._result.build()

    def add_dependencies(self, keywords):
        """Register input *DataStructure* objects as dependencies.

        Arguments:
            keywords (dict): User's keywords.
        """

        # no dependencies to add


CALC_MATR_ELEM = ComputeElementaryMatrix.run
