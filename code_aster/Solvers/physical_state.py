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
    """This object represents a Physical State of the model.

    Actually, it stores a stack of physical states et works as an *adapter*
    on the current state. Only the current state, the working one, has setters
    and so is writable. All other states on the stack are read-only.
    """

    class State:
        """Represents an elementary physical state (private)."""

        _time_prev = _time_step = None
        _primal_prev = _primal_step = _stress = _internVar = _externVar = None
        __setattr__ = no_new_attributes(object.__setattr__)

        @property
        def time_prev(self):
            """float: Previous time."""
            return self._time_prev

        @property
        def time_curr(self):
            """float: Current time."""
            return self._time_prev + self._time_step

        @property
        def time_step(self):
            """float: Time step."""
            return self._time_step

        @property
        def primal_prev(self):
            """FieldOnNodesReal: Primal field at previous time."""
            return self._primal_prev

        @property
        def primal_curr(self):
            """FieldOnNodesReal: Primal field at current time."""
            return self._primal_prev + self._primal_step

        @property
        def primal_step(self):
            """FieldOnNodesReal: Primal increment."""
            return self._primal_step

        @property
        def stress(self):
            """FieldOnCellsReal: Stress field."""
            return self._stress

        @property
        def internVar(self):
            """FieldOnCellsReal: Internal state variables."""
            return self._internVar

        @property
        def externVar(self):
            """FieldOnCellsReal: External state variables."""
            return self._externVar

        def copy(self, other):
            """Copy the content of an object into the current one.

            Arguments:
                other (PhysicalState.State): Object to be copied.

            Return:
                PhysicalState.State: Current object.
            """
            self._time_prev = other.time_prev
            self._time_step = other.time_step
            self._primal_prev = other.primal_prev and other.primal_prev.copy()
            self._primal_step = other.primal_step and other.primal_step.copy()
            self._stress = other.stress and other.stress.copy()
            self._internVar = other.internVar and other.internVar.copy()
            self._externVar = other.externVar and other.externVar.copy()
            return self

        def debugPrint(self, label=""):
            print(f"*** {label}Physical State at", self._time, flush=True)
            values = self._primal_prevv.getValues()
            print("* primal_prev ", sum(values) / len(values), flush=True)
            if self._primal_step:
                values = self._primal_step.getValues()
                print("* primal_step", sum(values) / len(values), flush=True)
            if self._stress:
                values = self._stress.getValues()
                print("* stress     ", sum(values) / len(values), flush=True)
            if self._internVar:
                values = self._internVar.getValues()
                print("* internVar  ", sum(values) / len(values), flush=True)

    provide = SOP.PhysicalState

    _current = _stack = _size = _stash = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, size=1):
        assert size > 0, f"invalid value ({size}) for 'size'"
        super().__init__()
        self._current = PhysicalState.State()
        self._stack = []
        self._size = size
        self._stash = None

    def getState(self, index=0):
        """Return a physical state by index (relative position).

        Arguments:
            index (int): 0 means the current one, -1 the previous one and so on.

        Returns:
            PhysicalState.State: physical state.
        """
        if index == 0:
            return self.current
        if index < -self._size or index > 0:
            raise IndexError(f"invalid value ({index}) for 'index': {-self._size} <= index <= 0")
        return self._stack[index]

    @property
    def current(self):
        """PhysicalState.State: The current physical state, the working one."""
        return self._current

    @property
    def time_prev(self):
        """float: Previous time."""
        return self.current.time_prev

    @time_prev.setter
    def time_prev(self, value):
        """Set previous time.

        Arguments:
            value (float): previous time
        """
        self.current._time_prev = value

    @property
    def time_curr(self):
        """float: Current time."""
        return self.current.time_prev + self.current.time_step

    @time_curr.setter
    def time_curr(self, value):
        """Set current time.

        Arguments:
            value (float): current time
        """
        self.current._time_step = value - self.current.time_prev

    @property
    def time_step(self):
        """float: Time step."""
        return self.current.time_step

    @time_step.setter
    def time_step(self, value):
        """Set time step.

        Arguments:
            value (float): time step
        """
        self.current._time_step = value

    @property
    def primal_prev(self):
        """FieldOnNodesReal: Primal field at previous time."""
        return self.current.primal_prev

    @primal_prev.setter
    def primal_prev(self, field):
        """Set previous primal field.

        Arguments:
           field (FieldOnNodesReal): primal
        """
        assert isinstance(field, FieldOnNodesReal), f"unexpected type: {field}"
        self.current._primal_prev = field

    @property
    def primal_curr(self):
        """FieldOnNodesReal: Primal field at current time."""
        return self.current._primal_prev + self.current.primal_step

    @primal_curr.setter
    def primal_curr(self, field):
        """Set primal increment.

        Arguments:
            field (FieldOnNodesReal): Primal increment
        """
        assert field is None or isinstance(field, FieldOnNodesReal), f"unexpected type: {field}"
        self.current._primal_step = field - self.current._primal_prev

    @property
    def primal_step(self):
        """FieldOnNodesReal: Primal increment."""
        return self.current.primal_step

    @primal_step.setter
    def primal_step(self, field):
        """Set primal increment.

        Arguments:
            field (FieldOnNodesReal): Primal increment
        """
        assert field is None or isinstance(field, FieldOnNodesReal), f"unexpected type: {field}"
        self.current._primal_step = field

    @property
    def stress(self):
        """FieldOnCellsReal: Stress field."""
        return self.current.stress

    @stress.setter
    def stress(self, field):
        """Set Stress field.

        Arguments:
            field (FieldOnCellsReal): Stress field
        """
        if field:
            assert isinstance(field, FieldOnCellsReal), f"unexpected type: {field}"
        self.current._stress = field

    @property
    def internVar(self):
        """FieldOnCellsReal: Internal state variables."""
        return self.current.internVar

    @internVar.setter
    def internVar(self, field):
        """Set Internal state variables.

        Arguments:
            field (FieldOnCellsReal): Internal state variables
        """
        if field:
            assert isinstance(field, FieldOnCellsReal), f"unexpected type: {field}"
        self.current._internVar = field

    @property
    def externVar(self):
        """FieldOnCellsReal: External state variables."""
        return self.current.externVar

    @externVar.setter
    def externVar(self, field):
        """Set external state variables.

        Arguments:
            field (FieldOnCellsReal): external state variables
        """
        assert field is None or isinstance(field, FieldOnCellsReal), f"unexpected type: {field}"
        self.current._externVar = field

    def stash(self):
        """Stores the object state to provide transactionality semantics."""
        self._stash = PhysicalState.State().copy(self.current)

    def revert(self):
        """Revert the object to its previous state."""
        assert self._stash, "stash is empty!"
        self.current.copy(self._stash)
        self._stash = None

    @profile
    def commit(self):
        """Commits the current changes and add the state on the stack."""
        # do not use '+=' to create a new object (and not modify previous values)
        current = self.current
        current._primal_prev = current._primal_prev
        if current._primal_step:
            current._primal_prev += current._primal_step
        current._primal_step.setValues(0.0)
        current._time_prev += current._time_step
        current._time_step = 0.0
        self._stash = None
        if len(self._stack) >= self._size:
            self._stack.pop(0)
        self._stack.append(current)
        self._current = PhysicalState.State().copy(current)

    @profile
    def getCurrentDelta(self):
        """Return the delta for the current state between it has been stashed.

        Returns:
            dict: Delta between states as returned by py:method:`as_dict`.
        """
        return self._states_difference(self._stash, self.current)

    @profile
    def getDeltaBetweenStates(self, index1, index2):
        """Return the delta between two states.

        The states are extracted using `getState(index)`.
        It returns the difference `getState(index2) - getState(index1)`.

        Arguments:
            index1 (int): Index of the first state.
            index2 (int): Index of the second state.

        Returns:
            dict: Delta between states as returned by py:method:`as_dict`.
        """
        return self._states_difference(self.getState(index1), self.getState(index2))

    @staticmethod
    def _states_difference(one, two):
        """Delta between states as returned by py:method:`as_dict`."""
        quantity, _ = two.primal_curr.getPhysicalQuantity().split("_")
        delta_primal = two.primal_curr - one.primal_curr

        ret = {"SIEF_ELGA": None, "VARI_ELGA": None, quantity: delta_primal}
        if one.stress and two.stress:
            ret["SIEF_ELGA"] = two.stress - one.stress
        if one.internVar and two.internVar:
            ret["VARI_ELGA"] = two.internVar - one.internVar

        return ret

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

        if phys_pb.isMechanical():
            type_field = "SIEF_R"
        else:
            type_field = "FLUX_R"

        return self.createFieldOnCells(phys_pb, "ELGA", type_field, value)

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
        self.time_prev = 0.0
        self.time_step = 0.0
        self.primal_prev = self.createPrimal(phys_pb, 0.0)
        self.primal_step = self.createPrimal(phys_pb, 0.0)
        if phys_pb.getBehaviourProperty():
            self.stress = self.createStress(phys_pb, 0.0)
            if phys_pb.isMechanical():
                self.internVar = self.createInternalVariablesNext(phys_pb, 0.0)
                self.externVar = None

    def as_dict(self):
        """Returns the fields as a dict.

        Returns:
            dict: Dict of fields.
        """
        current = self.current

        def getStoringName(field, default=" "):
            if field:
                quantity, _ = field.getPhysicalQuantity().split("_")
                if isinstance(field, (FieldOnNodesReal)):
                    loc = "NOEU"
                else:
                    loc = field.getLocalization()

                return quantity + "_" + loc

            return default

        return {
            getStoringName(current.stress): current.stress,
            getStoringName(current.internVar, "VARI_ELGA"): current.internVar,
            getStoringName(current.primal_prev).split("_")[0]: current.primal_prev,
        }

    def debugPrint(self, label=""):
        print(
            f"*** {label}Stack contains states for t =",
            [state.time_prev for state in self._stack],
            flush=True,
        )
