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

from ...Utilities import no_new_attributes, profile
from ...Objects import FieldOnNodesReal, FieldOnCellsReal, NonLinearResult


class PhysicalState:
    """This object represents a Physical State of the model."""

    _time = _time_step = None
    _displ = _displ_incr = _variP = _stress = None
    __setattr__ = no_new_attributes(object.__setattr__)

    @property
    def time(self):
        """float: Current time."""
        return self._time

    @time.setter
    def time(self, value):
        """Set time step.

        Arguments:
            value (float): time step
        """
        self._time = value

    @property
    def time_step(self):
        """float: Time step."""
        return self._time_step

    @time_step.setter
    def time_step(self, value):
        """Set time step.

        Arguments:
            value (float): time step
        """
        self._time_step = value

    @property
    def displ(self):
        """FieldOnNodes: Displacement field."""
        return self._displ

    @displ.setter
    def displ(self, field):
        """Set displacement field.

         Arguments:
            field (FieldOnNodesReal): displacement
        """
        self._displ = field

    @property
    def displ_incr(self):
        """FieldOnNodes: Displacement increment."""
        return self._displ_incr

    @displ_incr.setter
    def displ_incr(self, field):
        """Set displacement increment.

        Arguments:
            field (FieldOnNodes): Displacement increment
        """
        self._displ_incr = field

    @property
    def stress(self):
        """FieldOnCells: Stress field."""
        return self._stress

    @stress.setter
    def stress(self, field):
        """Set Stress field.

        Arguments:
            field (FieldOnCells): Stress field
        """
        self._stress = field

    @property
    def variP(self):
        """FieldOnNCells: Internal state variables."""
        return self._variP

    @variP.setter
    def variP(self, field):
        """Set Internal state variables.

        Arguments:
            field (FieldOnCells): Internal state variables
        """
        self._variP = field

    @profile
    def createDisplacement(self, phys_pb, value):
        """Create displacement field with a given value

        Arguments:
            phys_pb (PhysicalProblem): Physical problem
            value (float): value to set everywhere

        Returns:
            FieldOnNodes: displacement field with a given value (DEPL)
        """
        field = FieldOnNodesReal(phys_pb.getDOFNumbering())
        field.setValues(value)
        return field

    @profile
    def createFieldOnCells(self, phys_pb, type, value):
        """Create a field with a given value

        Arguments:
            phys_pb (PhysicalProblem): Physical problem
            type (str): type of the field to create
            value (float): value to set everywhere

        Returns:
            FieldOnCells: field with a given value of type "type"
        """
        field = FieldOnCellsReal(
            phys_pb.getModel(),
            phys_pb.getBehaviourProperty(),
            type,
            phys_pb.getElementaryCharacteristics(),
        )
        field.setValues(value)
        return field

    @profile
    def createStress(self, phys_pb, value):
        """Create stress field with a given value

        Arguments:
            phys_pb (PhysicalProblem): Physical problem
            value (float): value to set everywhere

        Returns:
            FieldOnCells: Stress field with a given value (SIEF_ELGA)
        """

        return self.createFieldOnCells(phys_pb, "ELGA_SIEF_R", value)

    @profile
    def createInternalVariablesNext(self, phys_pb, value):
        """Create internal state variables field with a given value

        Arguments:
            phys_pb (PhysicalProblem): Physical problem
            value (float): value to set everywhere

        Returns:
            FieldOnCells: internal state variables field with a given value (VARI_ELGA)
        """

        return self.createFieldOnCells(phys_pb, "ELGA_VARI_R", value)

    @profile
    def zeroInitialState(self, phys_pb):
        """Initialize with zero initial state

        Arguments:
            phys_pb (PhysicalProblem): Physical problem
        """
        self._displ = self.createDisplacement(phys_pb, 0.0)
        self._variP = self.createInternalVariablesNext(phys_pb, 0.0)
        self._stress = self.createStress(phys_pb, 0.0)
        self._time = 0.0
        self._time_step = 0.0

    @profile
    def readInitialState(self, phys_pb, params):
        """Initialize states

        Arguments:
            phys_pb (PhysicalProblem): Physical problem
            params (dict): dict of user's keywords
        """
        self.zeroInitialState(phys_pb)

        try:
            init_time = params.get("INCREMENT").get("LIST_INST").getValues()[0]
        except AttributeError:
            init_time = 0
        self._time = init_time

        init_params = params.get("ETAT_INIT")
        if init_params is not None:
            if "DEPL" in init_params:
                displ = init_params.get("DEPL")
                assert isinstance(displ, FieldOnNodesReal)
                self._displ = displ
            if "SIGM" in init_params:
                stress = init_params.get("SIGM")
                assert isinstance(stress, FieldOnCellsReal)
                self._stress = stress
            if "VARI" in init_params:
                variP = init_params.get("VARI")
                assert isinstance(variP, FieldOnCellsReal)
                self._variP = variP
            if "EVOL_NOLI" in init_params:
                resu = init_params.get("EVOL_NOLI")
                assert isinstance(resu, NonLinearResult)
                if "INST_ETAT_INIT" in init_params:
                    self._time = float(init_params.get("INST_ETAT_INIT"))
                    rank = resu.getNumberOfRanks() - 1
                    self.extractFieldsFromResult(resu, rank, ["DEPL", "SIEF_ELGA", "VARI_ELGA"])

    @profile
    def update(self, other):
        """Update current physical state with the previous one.

        Arguments:
            other (PhysicalState): physical model
        """
        # displ in two steps to create a new object (and not modified previous values)
        displ_up = self._displ + other.displ_incr
        self._displ = displ_up

        self._variP = other.variP
        self._stress = other.stress
        self._time += other.time_step

    @profile
    def extractFieldsFromResult(self, resu, rank, fields):
        """Extract fields at a rank

        Arguments:
            resu (NonLinearResult):
            rank (int): rank
            fields (str|tuple(str)): list of field to extract
        """

        if isinstance(fields, str):
            fields = tuple(fields)

        for field in fields:
            if field == "DEPL":
                self._displ = resu.getFieldOnNodesReal("DEPL", rank)
            elif field == "SIEF_ELGA":
                self._stress = resu.getFieldOnCellsReal("SIEF_ELGA", rank)
            elif field == "VARI_ELGA":
                self._variP = resu.getFieldOnCellsReal("VARI_ELGA", rank)
            else:
                raise RuntimeError("Unknown field")

    def as_dict(self):
        """Returns the fields as a dict.

        Returns:
            dict: Dict of fields.
        """
        return dict(DEPL=self._displ, SIEF_ELGA=self._stress, VARI_ELGA=self._variP)
