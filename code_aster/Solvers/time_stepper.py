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


class TimeStepper(SolverFeature):
    """This object deals with the time steps.

    Arguments:
        times (list[float]): List of time steps.
        epsilon (float): Value used to check equality between two times
            (default: 1.e-6).
    """

    provide = SOP.TimeStepper

    def __init__(self, times, epsilon=1e-6):
        super().__init__()
        times = list(times)
        assert sorted(times) == times
        self._times = times
        self._eps = epsilon
        self._current = 0
        self._first = 0
        self._last = len(times) - 1
        self._level = 0

    @property
    def null_increment(self):
        """float: Value of a null increment."""
        # FIXME could be just less than `self._eps`, no?
        return -1.0e150

    def setInitialStep(self, time):
        """Only use the times greater than `time`.
        If `time` is not already in the sequence, it is inserted.
        The initial time should be validated by `start()` to point to the next
        time step.

        Arguments:
            time (float): First time to be used.
        """
        times = self._times
        size = len(times)
        for idx in range(size):
            if time < times[idx] - self._eps:
                self._insert(idx, time)
                self._current = self._first = idx
                return
            elif abs(times[idx] - time) <= self._eps:
                self._current = self._first = idx
                return
        self._insert(size, time)
        self._current = self._first = size
        self._last = max(self._first, self._last)

    def setFinalStep(self, time=None):
        """Limit the sequence to the times lower than `time`.
        If `time` is not already in the sequence, it is inserted.
        If `time` is not provided, use the sequence up to its end.

        Arguments:
            time (float, optional): Last time to be used.
        """
        times = self._times
        time = times[-1] if time is None else time
        size = len(times)
        for idx in range(size):
            if time < times[idx] - self._eps:
                self._insert(idx, time)
                self._last = idx
                return
            elif abs(times[idx] - time) <= self._eps:
                self._last = idx
                return
        self._insert(size, time)
        self._last = size

    def size(self):
        """Return the number of steps to be calculated.

        Returns:
            int: Number of steps, including the initial one.
        """
        return self._last - self._first + 1

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
        if index <= self._first:
            self._first += 1
        if index <= self._current:
            # raise KeyError("can not insert a step before the current time.")
            self._current += 1
        self._times.insert(index, time)
        if index <= self._last:
            self._last += 1
        # print("insert at", index, self._first, self._last, self._current)

    def getInitialTime(self):
        """Returns the initial time.

        Returns:
            float: Initial time value.
        """
        return self._times[self._first]

    def getPrevious(self):
        """Returns the previous calculated step.

        Returns:
            float: Previous time value.
        """
        if self.hasFinished() or self._current - 1 < 0:
            raise IndexError

        return self._times[self._current - 1]

    def getCurrent(self):
        """Returns the current step, this to be calculated.

        Returns:
            float: Next time value.
        """
        if self.hasFinished():
            raise IndexError
        return self._times[self._current]

    def getIncrement(self):
        """Returns the increment to next step to be calculated.

        Returns:
            float: increment to the next time value.
        """
        if self._current == 0:
            return self.null_increment

        return self.getCurrent() - self.getPrevious()

    def start(self):
        """Start, move at the first step."""
        self._current = self._first + 1

    def completed(self):
        """Register the current step as completed successfully."""
        self._current += 1

    def split(self, nb_steps):
        """Split last time step in uniform sub-steps.

        Arguments:
            nb_steps (int): Number of sub-steps.
        """
        assert nb_steps > 1 and self._level < 20
        self._level += 1

        time_step = self.getIncrement() / nb_steps
        new = self.getCurrent()
        for _ in range(max(0, nb_steps - 1)):
            new -= time_step
            self._insert(self._current, new)
            self._current -= 1

    def raiseError(self, exc):
        """Raise an error executing the last step.

        Arguments:
            exc (Exception): Error to be raised.
        """
        # manage substepping...
        # for example, step back to the previous step
        self._current = max(self._current - 1, 0)
        raise exc

    def _debug_(self, prefix=""):
        prefix = prefix or "TimeStepper"
        print(prefix, self._first, self._last, self._current, self._times)
