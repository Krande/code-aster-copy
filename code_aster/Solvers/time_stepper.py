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

from .solver_features import SolverFeature
from .solver_features import SolverOptions as SOP
from ..Utilities import no_new_attributes, logger


class TimeStepper(SolverFeature):
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

    _times = _eps = _current = _initial = _final = _last = _level = None
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
        self._level = 0

    @property
    def null_increment(self):
        """float: delta between two steps to be considered as null."""
        return self._eps

    def _check_bounds(self):
        """Remove out of bounds values."""
        times = self._times
        if self._initial is not None:
            while times and times[0] < self._initial + self._eps:
                times.pop(0)
        if self._final is not None:
            while times and times[-1] > self._final + self._eps:
                times.pop()
            if not times or self._final > times[-1] + self._eps:
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

        The next calculated time is reset to the first in the list.
        Calling `setInitial` resets the current step at the first position.

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
        Calling `setFinal` resets the current step at the first position.

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

    def split(self, nb_steps):
        """Split last time step in uniform sub-steps.

        Arguments:
            nb_steps (int): Number of sub-steps.
        """
        assert nb_steps > 1 and self._level < 20
        # not a splitting level, but a number of splits
        self._level += 1

        time_step = self.getIncrement() / nb_steps
        new = self.getCurrent()
        for _ in range(max(0, nb_steps - 1)):
            new -= time_step
            self._insert(self._current, new)

    def raiseError(self, exc):
        """Raise an error executing the last step.

        Arguments:
            exc (Exception): Error to be raised.
        """
        # manage substepping...
        # for example, step back to the previous step
        self._current = max(self._current - 1, 0)
        raise exc

    def __repr__(self):
        return f"<TimeStepper(from {self._initial} to {self._final}, size {self.size()}: {self._times})>"

    @staticmethod
    def from_keywords(**args):
        """Initialize a TimeStepper from user keywords as provided to
        the INCREMENT common set of keywords.

        Arguments:
            args (dict): keywords as for INCREMENT.

        Returns:
            TimeStepper: a new TimeStepper.
        """
        assert "LIST_INST" in args, "THER_NON_LINE not yet supported!"
        times = args["LIST_INST"].getValues()
        initial = times[0]
        if "INST_INIT" in args:  # because None has a special meaning
            initial = args["INST_INIT"]
        if args.get("NUME_INST_INIT"):
            initial = times[args["NUME_INST_INIT"]]
        final = args.get("INST_FIN")
        if args.get("NUME_INST_FIN"):
            final = times[args["NUME_INST_FIN"]]
        stp = TimeStepper(times, initial=initial, final=final, epsilon=args["PRECISION"])
        return stp
