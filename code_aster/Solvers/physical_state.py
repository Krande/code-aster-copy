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

from ..Objects import DiscreteComputation, FieldOnCellsReal, FieldOnNodesReal
from ..Utilities import no_new_attributes, profile
from .base_features import BaseFeature
from .solver_features import SolverOptions as SOP


class PhysicalState(BaseFeature):
    """This object represents a Physical State of the model."""

    provide = SOP.PhysicalState

    _time = _time_step = None
    _primal = _primal_step = _internVar = _stress = None
    _externVar = _externVar_next = None
    _stash = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        super().__init__()
        self._stash = None

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

    def copy(self, other):
        """Copy the content of an object into the current one.

        Arguments:
            other (PhysicalState): Object to be copied.
        """
        self._time = other.time
        self._time_step = other.time_step
        self._primal = other.primal and other.primal.duplicate()
        self._primal_step = other.primal_step and other.primal_step.duplicate()
        self._stress = other.stress and other.stress.duplicate()
        self._internVar = other.internVar and other.internVar.duplicate()
        self._externVar = other.externVar and other.externVar.duplicate()
        self._externVar_next = other.externVar_next and other.externVar_next.duplicate()

    def stash(self):
        """Stores the object state to provide transactionality semantics."""
        self._stash = PhysicalState()
        self._stash.copy(self)

    @profile
    def getIncrement(self):
        """Return the delta between the previous and the current state.

        The previous state is expected to be found in the stash.

        Returns:
            dict: Delta between states as returned by py:method:`as_dict`.
        """
        quantity, _ = self._primal.getPhysicalQuantity().split("_")
        return {
            "SIEF_ELGA": self._stress - self._stash._stress,
            "VARI_ELGA": self.internVar - self._stash._internVar,
            quantity: self._primal + self._primal_step - self._stash._primal,
        }

    def revert(self):
        """Revert the object to its previous state."""
        assert self._stash, "stash is empty!"
        self.copy(self._stash)
        self._stash = None

    @profile
    def commit(self):
        """Commits the current changes."""
        # do not use '+=' to create a new object (and not modified previous values)
        self._primal = self._primal + self._primal_step
        self._primal_step = None
        self._time += self._time_step
        self._time_step = 0.0
        self._stash = None

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
        self._time = 0.0
        self._time_step = 0.0
        self._primal = self.createPrimal(phys_pb, 0.0)
        self._primal_step = None
        self._stress = self.createStress(phys_pb, 0.0)
        self._internVar = self.createInternalVariablesNext(phys_pb, 0.0)
        self._externVar = None
        self._externVar_next = None

    def as_dict(self):
        """Returns the fields as a dict.

        Returns:
            dict: Dict of fields.
        """
        quantity, _ = self._primal.getPhysicalQuantity().split("_")
        return {"SIEF_ELGA": self._stress, "VARI_ELGA": self.internVar, quantity: self._primal}

    def debugPrint(self, label=""):
        print(f"*** {label}Physical State at", self._time, flush=True)
        values = self._primal.getValues()
        print("* primal     ", sum(values) / len(values), flush=True)
        if self._primal_step:
            values = self._primal_step.getValues()
            print("* primal_step", sum(values) / len(values), flush=True)
        if self._stress:
            values = self._stress.getValues()
            print("* stress     ", sum(values) / len(values), flush=True)
        if self._internVar:
            values = self._internVar.getValues()
            print("* internVar  ", sum(values) / len(values), flush=True)
