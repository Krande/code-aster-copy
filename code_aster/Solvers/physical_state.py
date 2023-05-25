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

from ..Objects import DiscreteComputation, FieldOnCellsReal, FieldOnNodesReal, NonLinearResult
from ..Utilities import no_new_attributes, profile
from .base_features import BaseFeature
from .solver_features import SolverOptions as SOP
from ..Commands import IMPR_CO


class PhysicalState(BaseFeature):
    """This object represents a Physical State of the model."""

    provide = SOP.PhysicalState

    _time = _time_step = None
    _primal = _primal_step = _internVar = _stress = None
    _externVar = _externVar_next = None
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
    def primal(self):
        """FieldOnNodesReal: Primal field."""
        return self._primal

    @primal.setter
    def primal(self, field):
        """Set primal field.

        Arguments:
           field (FieldOnNodesReal): primal
        """
        assert isinstance(field, FieldOnNodesReal), f"unexpected type: {field}"
        self._primal = field

    @property
    def primal_step(self):
        """FieldOnNodesReal: Primal increment."""
        return self._primal_step

    @primal_step.setter
    def primal_step(self, field):
        """Set primal increment.

        Arguments:
            field (FieldOnNodesReal): Primal increment
        """
        self._primal_step = field

    @property
    def stress(self):
        """FieldOnCellsReal: Stress field."""
        return self._stress

    @stress.setter
    def stress(self, field):
        """Set Stress field.

        Arguments:
            field (FieldOnCellsReal): Stress field
        """
        assert isinstance(field, FieldOnCellsReal), f"unexpected type: {field}"
        self._stress = field

    @property
    def internVar(self):
        """FieldOnCellsReal: Internal state variables."""
        return self._internVar

    @internVar.setter
    def internVar(self, field):
        """Set Internal state variables.

        Arguments:
            field (FieldOnCellsReal): Internal state variables
        """
        assert isinstance(field, FieldOnCellsReal), f"unexpected type: {field}"
        self._internVar = field

    @property
    def externVar(self):
        """FieldOnCellsReal: External state variables."""
        return self._externVar

    @externVar.setter
    def externVar(self, field):
        """Set external state variables.

        Arguments:
            field (FieldOnCellsReal): external state variables
        """
        assert isinstance(field, FieldOnCellsReal), f"unexpected type: {field}"
        self._externVar = field

    @property
    def externVar_next(self):
        """FieldOnCellsReal: External state variables at end of step."""
        return self._externVar_next

    @externVar_next.setter
    def externVar_next(self, field):
        """Set external state variables at end of step

        Arguments:
            field (FieldOnCellsReal): external state variables at end of step
        """
        self._externVar_next = field

    # FIXME setPrimalValue?
    @profile
    def createPrimal(self, phys_pb, value):
        """Create primal field with a given value

        Arguments:
            phys_pb (PhysicalProblem): Physical problem
            value (float): value to set everywhere

        Returns:
            FieldOnNodes: primal field with a given value (DEPL|TEMP)
        """
        field = FieldOnNodesReal(phys_pb.getDOFNumbering())
        field.setValues(value)
        return field

    @profile
    def createFieldOnCells(self, phys_pb, localization, quantity, value):
        """Create a field with a given value

        Arguments:
            phys_pb (PhysicalProblem): Physical problem
            localization (str): localization of the field to create (ELNO, ELEM, ELGA)
            quantity (str): type of the field to create (ex: SIEF_R, ...)
            value (float): value to set everywhere

        Returns:
            FieldOnCells: field with a given value of type "type"
        """
        assert phys_pb.getBehaviourProperty(), "unexpected empty BehaviourProperty"
        field = FieldOnCellsReal(
            phys_pb.getModel(),
            localization,
            quantity,
            phys_pb.getBehaviourProperty(),
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
        return self.createFieldOnCells(phys_pb, "ELGA", "SIEF_R", value)

    @profile
    def createInternalVariablesNext(self, phys_pb, value):
        """Create internal state variables field with a given value

        Arguments:
            phys_pb (PhysicalProblem): Physical problem
            value (float): value to set everywhere

        Returns:
            FieldOnCells: internal state variables field with a given value (VARI_ELGA)
        """
        return self.createFieldOnCells(phys_pb, "ELGA", "VARI_R", value)

    @profile
    def createTimeField(self, phys_pb, value):
        """Create time field with a given value

        Arguments:
            phys_pb (PhysicalProblem): Physical problem
            value (float): value to set everywhere

        Returns:
            ConstantFieldOnCellsReal: time field with a given value
        """
        disc_comp = DiscreteComputation(phys_pb)
        return disc_comp.createTimeField(value)

    @profile
    def zeroInitialState(self, phys_pb):
        """Initialize with zero initial state

        Arguments:
            phys_pb (PhysicalProblem): Physical problem
        """
        self._primal = self.createPrimal(phys_pb, 0.0)
        self._internVar = self.createInternalVariablesNext(phys_pb, 0.0)
        self._stress = self.createStress(phys_pb, 0.0)
        self._time = 0.0
        self._time_step = 0.0

    # FIXME set 'other' optional? removed?
    @profile
    def update(self, other):
        """Update current physical state with the previous one.

        Arguments:
            other (PhysicalState): physical model
        """
        # primal in two steps to create a new object (and not modified previous values)
        primal_up = self._primal + other.primal_step
        self._primal = primal_up
        self._internVar = other.internVar
        self._externVar = other.externVar_next
        self._stress = other.stress
        self._time += other.time_step

    def as_dict(self):
        """Returns the fields as a dict.

        Returns:
            dict: Dict of fields.
        """
        quantity, fld_type = self._primal.getPhysicalQuantity().split("_")
        return {"SIEF_ELGA": self._stress, "VARI_ELGA": self.internVar, quantity: self._primal}

    def debugPrint(self):
        if len(self._primal.getName()) > 8:
            IMPR_CO(CHAINE=self._primal.getName() + ".VALE", NIVEAU=-1, UNITE=6)
        else:
            IMPR_CO(CONCEPT=_F(NOM=self._primal), NIVEAU=-1, UNITE=6)
        if self._primal_step:
            if len(self._primal_step.getName()) > 8:
                IMPR_CO(CHAINE=self._primal_step.getName() + ".VALE", NIVEAU=-1, UNITE=6)
            else:
                IMPR_CO(CONCEPT=_F(NOM=self._primal_step), NIVEAU=-1, UNITE=6)
        if self._stress:
            if len(self._stress.getName()) > 8:
                IMPR_CO(CHAINE=self._stress.getName() + ".CELV", NIVEAU=-1, UNITE=6)
            else:
                IMPR_CO(CONCEPT=_F(NOM=self._stress), NIVEAU=-1, UNITE=6)
        if self._internVar:
            if len(self._internVar.getName()) > 8:
                IMPR_CO(CHAINE=self._internVar.getName() + ".CELV", NIVEAU=-1, UNITE=6)
            else:
                IMPR_CO(CONCEPT=_F(NOM=self._internVar), NIVEAU=-1, UNITE=6)
