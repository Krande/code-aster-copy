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

import unittest

from code_aster.Solvers.nl_integr_features import NLIntegrFeature
from code_aster.Solvers.nl_integr_features import NLIntegrOptions as SOP


class TimeStepper(NLIntegrFeature):
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


class BasicTest(unittest.TestCase):
    """Check for internal methods."""

    def test_initial(self):
        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        stp.setInitialStep(0.0)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.0)
        stp.start()
        self.assertEqual(stp.size(), 4)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.25)

        eps = 1.0e-3
        stp = TimeStepper([0.25, 1.0, 2.0], eps)
        stp.setInitialStep(0.0)
        self.assertEqual(stp.size(), 4)
        stp.start()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.25)

        stp.setInitialStep(0.3)
        self.assertEqual(stp.size(), 3)
        stp.start()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)

        stp.setInitialStep(0.1)
        stp.start()
        self.assertEqual(stp.size(), 5)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.25)

        stp.setInitialStep(0.1 + eps * 0.999)
        stp.start()
        self.assertEqual(stp.size(), 5)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.25)

        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        stp.setInitialStep(1.0)
        self.assertEqual(stp.size(), 2)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)

        stp.setInitialStep(2.5)
        self.assertEqual(stp.size(), 1)
        self.assertFalse(stp.hasFinished())
        stp.start()
        self.assertTrue(stp.hasFinished())

    def test_final(self):
        eps = 1.0e-3
        stp = TimeStepper([0.0, 0.25, 1.0, 2.0], eps)
        stp.setFinalStep(2.0 - eps * 0.999)
        self.assertEqual(stp.size(), 4)
        stp.setFinalStep(2.0 + eps * 0.999)
        self.assertEqual(stp.size(), 4)

        stp.setFinalStep(1.9)
        self.assertEqual(stp.size(), 4)
        for _ in range(3):
            stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.9)

        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        stp.setFinalStep(2.5)
        self.assertEqual(stp.size(), 5)
        for _ in range(4):
            stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.5)

        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        self.assertEqual(stp.size(), 4)
        stp.setFinalStep(0.8)
        self.assertEqual(stp.size(), 3)
        stp.setFinalStep()
        self.assertEqual(stp.size(), 5)

    def test_basic(self):
        stp = TimeStepper([0.0, 0.25, 1.0, 2.0, 3.0])
        self.assertFalse(stp.hasFinished())

        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.0)
        delta_t = stp.getIncrement()
        self.assertEqual(delta_t, stp.null_increment)

        stp.setInitialStep(0.0)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.0)
        stp.start()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.25)
        delta_t = stp.getIncrement()
        self.assertAlmostEqual(delta_t, 0.25)
        stp.completed()

        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        delta_t = stp.getIncrement()
        self.assertAlmostEqual(delta_t, 0.75)
        stp.completed()

        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)

        with self.assertRaises(ValueError):
            stp.raiseError(ValueError())

        self.assertFalse(stp.hasFinished())
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)
        stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 3.0)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 3.0)
        self.assertFalse(stp.hasFinished())

    def test_restart(self):
        stp = TimeStepper([0.0, 0.25, 1.0, 2.0, 3.0])
        self.assertEqual(stp.size(), 5)
        stp.setFinalStep(1.0)
        self.assertEqual(stp.size(), 3)

        stp.start()
        self.assertAlmostEqual(stp.getCurrent(), 0.25)
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 1.0)
        stp.completed()
        self.assertTrue(stp.hasFinished())

        stp.setInitialStep(1.0)
        self.assertEqual(stp.size(), 1)
        stp.setFinalStep(3.0)
        self.assertEqual(stp.size(), 3)
        stp.start()
        self.assertAlmostEqual(stp.getCurrent(), 2.0)

    def test_split(self):
        stp = TimeStepper([0.0, 1.0, 2.0, 3.0])
        self.assertEqual(stp.size(), 4)
        stp.start()
        self.assertAlmostEqual(stp.getCurrent(), 1.0)
        stp.split(2)
        self.assertEqual(stp.size(), 5)
        self.assertAlmostEqual(stp.getCurrent(), 0.5)
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 1.0)
        stp.split(5)
        self.assertEqual(stp.size(), 9)
        self.assertAlmostEqual(stp.getCurrent(), 0.6)
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 0.7)


if __name__ == "__main__":
    unittest.main()
