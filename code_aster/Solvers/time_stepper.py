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

from enum import Flag, auto
from math import log10

import numpy
from libaster import ConvergenceError, IntegrationError, SolverError

from ..Cata.Syntax import _F
from ..Messages import MessageLog
from ..Utilities import force_list, logger, no_new_attributes
from .base_features import Observer
from .solver_features import SolverFeature
from .solver_features import SolverOptions as SOP


class TimeStepper(SolverFeature, Observer):
    """This object deals with the time steps.

    It gives the list of the time steps to be calculated. The initial time is
    not in the list since it is already known. The final time is.
    The initial time may be *None*, undeterminated. In this case, the first
    increment is also *None*.

    Arguments:
        times (list[float]): List of time steps to be calculated.
        epsilon (float, optional): Value used to check equality between two times
            (default: 1.e-6).
        initial (float, optional): Initial time (default: 0.0).
        final (float, optional): Final time (default: the last given).
    """

    provide = SOP.TimeStepper

    _times = _eps = _current = _initial = _final = _last = None
    _actions = _state = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, times, epsilon=1e-6, initial=0.0, final=None):
        super().__init__()
        times = list(times)
        if sorted(times) != times:
            raise ValueError("the time steps must be ordered")
        self._times = times
        self._eps = epsilon
        self._initial = initial
        self._final = final
        logger.debug("TimeStepper.init: %s, %s, %s", initial, final, times)
        self._check_bounds()
        self._actions = []
        # FIXME the previous phys_state should be hold by NonLinearSolver
        self._state = {}

    @property
    def null_increment(self):
        """float: delta between two steps to be considered as null."""
        return self._eps

    def copy(self):
        """Return a copy of the object.

        Returns:
            TimeStepper: copy of the object.
        """
        new = TimeStepper(self._times, initial=self._initial, final=self._final, epsilon=self._eps)
        for act in self._actions:
            new.register_event(act.copy())
        new.register_default_error_event()
        return new

    def _check_bounds(self):
        """Remove out of bounds values."""
        times = self._times
        if self._initial is not None:
            while times and self.cmp(times[0], self._initial) <= 0:
                times.pop(0)
        if self._final is not None:
            while times and self.cmp(times[-1], self._final) > 0:
                times.pop()
            if not times or self.cmp(self._final, times[-1]) > 0:
                times.append(self._final)
        self._current = 0
        if not times:
            # empty list
            self._last = -1
            self._final = self._initial
            return
        self._last = len(times) - 1
        self._final = times[-1]

    def setInitial(self, time):
        """Define the initial time. Lesser values are removed.

        Calling `setInitial` resets the next calculated time at the first
        position.

        Arguments:
            time (float): First time to be used.
        """
        self._initial = time
        self._check_bounds()

    def setFinal(self, time=None):
        """Limit the sequence to the times lower than `time`.

        If `time` is not already in the sequence, it is appended.
        If `time` is not provided, the final time is set to the last one of
        the sequence.
        Calling `setInitial` resets the next calculated time at the first
        position.

        Arguments:
            time (float, optional): Last time to be used.
        """
        self._final = time
        self._check_bounds()

    def size(self):
        """Return the total number of steps in the list.

        Returns:
            int: Number of steps.
        """
        return len(self._times)

    def remaining(self):
        """Return the number of steps not yet calculated.

        Returns:
            int: Number of steps.
        """
        return self._last - self._current + 1

    def hasFinished(self):
        """Tell if there are steps to be computed.

        Returns:
            bool: *True* if there is no step to be computed, *False* otherwise.
        """
        return self._current > self._last

    def _insert(self, index, time):
        """Inserts the step to given index. The caller must check for already
        existing time.

        Arguments:
            index (int): index to insert nex step
            time (float): time value to insert.
        """
        # print("\ninsert at", index, time, self._current, end=" ")
        if index < self._current:
            raise KeyError("can not insert a step before the current time.")
        self._times.insert(index, time)
        self._last = len(self._times) - 1
        # print("\n->", self._current, self._last, self._times, flush=True)

    def getInitial(self):
        """Returns the initial time (not calculated).

        Returns:
            float: Initial time value.
        """
        return self._initial

    def getFinal(self):
        """Returns the last time to be calculated.

        Returns:
            float: Final time value.
        """
        return self._final

    def getPrevious(self):
        """Returns the previous calculated step.

        Returns:
            float: Previous time value.
        """
        if self.hasFinished():
            raise IndexError("no more timesteps")
        if self._current == 0:
            return self._initial

        return self._times[self._current - 1]

    def getCurrent(self):
        """Returns the current step, this to be calculated.

        Returns:
            float: Next time value.
        """
        if self.hasFinished():
            raise IndexError("no more timesteps")
        return self._times[self._current]

    def getIncrement(self):
        """Returns the increment to next step to be calculated.

        Returns:
            float: increment to the next time value.
        """
        prev = self.getInitial() if self._current == 0 else self.getPrevious()
        if prev is None:
            return None
        return self.getCurrent() - prev

    def completed(self):
        """Register the current step as completed successfully."""
        if self.hasFinished():
            raise IndexError("no more timesteps")
        self._current += 1

    def update(self, event):
        """Receive update from event.

        Arguments:
            event (EventSource): Object that sends the notification.
        """
        eid, data = event.get_state()
        if eid & SOP.IncrementalSolver:
            logger.info("#DBG from incr: %s", data)

    def cmp(self, time1, time2):
        """Compare two times using epsilon.

        Arguments:
            time1 (float): first argument.
            time2 (float): second argument.

        Returns:
            int: -1 if time1 < time2, 0 if time1 == time2, +1 if time1 > time2
            using epsilon.
        """
        return (time1 > time2 + self._eps) - (time1 + self._eps < time2)

    def __repr__(self):
        return f"<TimeStepper(from {self._initial} to {self._final}, size {self.size()}: {self._times})>"

    @staticmethod
    def from_keywords(**args):
        """Initialize a TimeStepper from user keywords as provided to
        the INCREMENT common set of keywords.

        Arguments:
            args (dict): keywords as for INCREMENT.

        Returns:
            TimeStepper: a new TimeStepper or a copy of the object from INCREMENT.
        """
        assert "LIST_INST" in args, "THER_NON_LINE not yet supported!"
        try:
            stp = args["LIST_INST"].stepper.copy()
            times = stp._times
        except AttributeError:
            # ListOfFloats
            times = args["LIST_INST"].getValues()
            stp = None
        initial = times[0]
        eps = args.get("PRECISION", 1.0e-6)
        if "INST_INIT" in args:  # because None has a special meaning
            initial = args["INST_INIT"]
        if args.get("NUME_INST_INIT"):
            initial = times[args["NUME_INST_INIT"]]
        final = args.get("INST_FIN")
        if args.get("NUME_INST_FIN"):
            final = times[args["NUME_INST_FIN"]]
        if stp is None:
            stp = TimeStepper(times, initial=initial, final=final, epsilon=eps)
        else:
            stp.setInitial(initial)
            stp.setFinal(final)
            stp._eps = eps
        stp.register_default_error_event()
        return stp

    @staticmethod
    def command_factory(args):
        """Create a TimeStepper from DEFI_LIST_INST keywords.

        *Transitional function during migration from legacy operators that need
        a TimesList object and the ones are used a TimeStepper.*

        Argumentss:
            args (dict): User keywords

        Returns:
            TimeStepper: a new TimeStepper.
        """
        args = _F(args)
        definition = args["DEFI_LIST"]
        if "VALE" in definition:
            times = definition["VALE"]
        elif "LIST_INST" in definition:
            times = definition["LIST_INST"].getValues()
        else:
            # this option did not exist with AUTO
            result = definition["RESULTAT"]
            div = definition["SUBD_PAS"]
            orig = [result.getTime(idx) for idx in result.getIndexes()]
            times = [orig.pop(0)]
            for step in orig:
                times.extend(numpy.linspace(times[-1], step, div + 1)[1:])
        stp = TimeStepper(times, initial=None)
        maxLevel = 0
        for fail in args["ECHEC"]:
            if fail["EVENEMENT"] == "ERREUR":
                event = TimeStepper.Error()
            elif fail["EVENEMENT"] == "DELTA_GRANDEUR":
                event = TimeStepper.MaximumIncrement(
                    fieldName=fail["NOM_CHAM"],
                    component=fail["NOM_CMP"],
                    maxValue=fail["VALE_REF"],
                    group=None,
                )
            if not event:  # not yet supported, ignored
                continue
            if fail["ACTION"] == "ARRET":
                act = TimeStepper.Interrupt(event)
            elif fail["ACTION"] == "DECOUPE":
                if fail["SUBD_METHODE"] == "MANUEL":
                    maxLevel = max(maxLevel, fail["SUBD_NIVEAU"])
                    act = TimeStepper.Split(
                        event,
                        nbSteps=fail["SUBD_PAS"],
                        maxLevel=maxLevel,
                        minStep=fail["SUBD_PAS_MINI"],
                    )
                else:
                    assert fail["SUBD_METHODE"] == "AUTO"
                    act = TimeStepper.Split(event, minStep=fail["SUBD_PAS_MINI"])
            else:  # not yet supported, ignored
                continue
            stp.register_event(act)
        stp.register_default_error_event()
        for act in stp._actions:
            try:
                act.setMaxLevel(maxLevel)
            except AttributeError:
                pass
        return stp

    # event management
    def register_event(self, action):
        """Register an action to react to an event.

        Args:
            action (~TimeStepper.Action): Type of action.
        """
        self._actions.append(action)

    def register_default_error_event(self):
        """Register a default action for Error event if there is no one."""
        if not [act for act in self._actions if isinstance(act.event, TimeStepper.Error)]:
            self.register_event(TimeStepper.Split(TimeStepper.Error(), nbSteps=2, maxLevel=4))

    def failed(self, exc):
        """React to a raised exception.

        Args:
            exc (AsterError): The exception just raised.
        """
        for act in self._actions:
            if not isinstance(act.event, TimeStepper.Error):
                continue
            if act.event.is_raised(exception=exc):
                act.call(timeStepper=self, exception=exc)
                return
        raise TypeError("should not pass here!")  # because of default error

    def check_event(self, phys_state):
        """Check if an event should be raised.

        Arguments:
            phys_state (PhysicalState): Current physical state.

        Returns:
            bool: *False* if something went wrong, *True** if the step is ok.
        """
        # compute increment
        prev = self._state.get("previous")
        delta = phys_state.as_dict()
        if prev:
            for name, field in delta.items():
                delta[name] - prev[name]
        self._state["previous"] = phys_state.as_dict()

        for act in self._actions:
            if not isinstance(act.event, TimeStepper.MaximumIncrement):
                continue
            if act.event.is_raised(delta=delta):
                act.call(timeStepper=self)
                # the exit code of 'call()' could be be ok (ADAPTATION?)...
                return False
        return True

    class Event:
        """Event raised in case of error."""

        # TODO
        # ResidueDivergence: raised if resi_t+2 > min(resi_t, resi_t+1) (DIVE_RESI)
        # MaximumResidue: raised if RESI_GLOB_MAXI > value before the end of iterations (RESI_MAXI)
        # COLLISION, INTERPENETRATION, INSTABILITE

        def is_raised(self, **context):
            """Tell if the event is raised in this context.

            Returns:
                bool: *True* if the event is raised, *False* otherwise.
            """
            raise NotImplementedError("must be subclassed!")

    class Error(Event):
        """Event raised in case of generic error."""

        @staticmethod
        def is_raised(**context):
            """Tell if the event is raised in this context.

            Returns:
                bool: *True* if the event is raised, *False* otherwise.
            """
            return isinstance(context.get("exception"), (ConvergenceError, IntegrationError))

    class MaximumIncrement(Event):
        """Event raised when the increment of a component exceeds a value
        (DELTA_GRANDEUR keyword).

        Arguments:
            fieldName (str): Name of the field to be checked.
            component (str): Component name.
            maxValue (float): Maximum increment.
            group (*misc*): Restrict the checking to a part of the model.
        """

        def __init__(self, fieldName, component, maxValue, group=None):
            self._fieldName = fieldName
            self._cmp = component
            self._value = maxValue
            self._group = force_list(group or [])

        def is_raised(self, **context):
            """Tell if the event is raised in this context.

            Returns:
                bool: *True* if the event is raised, *False* otherwise.
            """
            delta = context.get("delta")
            if not delta:
                return False
            field = delta.get(self._fieldName)
            if not field:
                return False
            maxIncr = max(field.getValuesWithDescription(self._cmp, self._group)[0])
            logger.info("#DBG check delta %s >? %s", maxIncr, self._value)
            raised = maxIncr > self._value
            if raised:
                logger.info(MessageLog.GetText("I", "MECANONLINE10_24"))
            return raised

    class Action:
        """This object allow to execute an action on a TimeStepper when an event
        occurs."""

        _event = None
        __setattr__ = no_new_attributes(object.__setattr__)

        def __init__(self, event) -> None:
            self._event = event

        def copy(self):
            """Return a copy of the object."""
            return self.__class__(self._event)

        @property
        def event(self):
            """TimeStepper.Event: event to react."""
            return self._event

        def call(self, **context):
            """Execute the action.

            Arguments:
                context (dict): Context of the event.
            """
            raise NotImplementedError("must be subclassed!")

    class Interrupt(Action):
        """This action stops the calculation (keyword value: ARRET)."""

        def call(self, **context):
            """Execute the action.

            Arguments:
                context (dict): Context of the event.
            """
            logger.info(MessageLog.GetText("I", "MECANONLINE10_30"))
            raise context["exception"]

    class Split(Action):
        """This action adds intermediate timesteps (keyword value: DECOUPE).

        Arguments:
            nbSteps (int): Number of sub-steps.
            maxLevel (int): Maximal number of levels.
            minStep (float): Minimal value of a step.
        """

        _nbSteps = _maxLevel = _minStep = None
        _reset = _stop = None

        def __init__(self, event, nbSteps, maxLevel=1e6, minStep=1.0e-12):
            super().__init__(event)
            assert nbSteps > 1, nbSteps
            self._nbSteps = nbSteps
            self._maxLevel = maxLevel
            self._minStep = minStep
            self._reset = []
            self._stop = TimeStepper.Interrupt(event)

        def copy(self):
            """Return a copy of the object."""
            return self.__class__(self._event, self._nbSteps, self._maxLevel, self._minStep)

        def setMaxLevel(self, maxLevel):
            """Set the maximum number of levels.

            Arguments:
                maxLevel (int): Maximal number of levels.
            """
            self._maxLevel = maxLevel

        def call(self, **context):
            """Execute the action.

            Arguments:
                context (dict): Context of the event.
            """
            logger.info(MessageLog.GetText("I", "MECANONLINE10_33"))
            logger.info(MessageLog.GetText("I", "SUBDIVISE_1"))
            stp = context["timeStepper"]
            last = stp.getCurrent()
            # TODO check that minStep > epsilon
            logger.debug("check splitting level: %s %s", self._reset, last)
            while self._reset and stp.cmp(last, self._reset[-1]) >= 0:
                self._reset.pop()
            self._reset.append(last)
            # not a splitting level, but a number of splits
            if len(self._reset) > self._maxLevel:
                self._stop.call(exception=SolverError("SUBDIVISE_17", (), (self._maxLevel,), ()))

            time_step = stp.getIncrement() / self._nbSteps
            if time_step < self._minStep:
                self._stop.call(exception=SolverError("SUBDIVISE_53", (), (), (self._minStep,)))
            logger.info(MessageLog.GetText("I", "SUBDIVISE_51"))
            new = stp.getCurrent()
            for _ in range(max(0, self._nbSteps - 1)):
                new -= time_step
                stp._insert(stp._current, new)

    class AutoSplit(Action):
        """This action adds intermediate timesteps (keyword value: DECOUPE)
        using less parameters than the *Split* action.

        Arguments:
            minStep (float): Minimal value of a step.
        """

        _minStep = None
        _manual = None

        def __init__(self, event, minStep):
            super().__init__(event)
            self._minStep = minStep
            self._manual = TimeStepper.Split(event, nbSteps=4, minStep=minStep)

        def copy(self):
            """Return a copy of the object."""
            return self.__class__(self._event, self._minStep)

        def call(self, **context):
            """Execute the action.

            Arguments:
                context (dict): Context of the event.
            """
            logger.info(MessageLog.GetText("I", "MECANONLINE10_33"))
            logger.info(MessageLog.GetText("I", "SUBDIVISE_2"))
            stp = context["timeStepper"]
            residuals = context["residuals"]
            xn = 0.0
            sx = 0.0
            sy = 0.0
            sxx = 0.0
            syx = 0.0
            for i, resi in enumerate(residuals):
                xx = log10(resi)
                if i >= len(residuals) - 2:
                    weight = 2.0
                else:
                    weight = 1.0
                xn += weight
                sx += weight * xx
                sy += weight * (i + 1)
                sxx += weight * xx * 2
                syx += weight * xx * (i + 1)
            a1 = sxx * sy - sx * syx
            a2 = -sx * sy + syx * xn
            b = -(sx**2) + sxx * xn
            assert False, "WTF"

    # ITER_SUPPL
    # AUTRE_PILOTAGE
    # ADAPT_COEF_PENA
    # CONTINUE
