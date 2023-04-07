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

from ...NonLinear import SolverFeature
from ...NonLinear import SolverOptions as SOP


class TimeStepper(SolverFeature):
    """ "Basic time stepper.

    Arguments:
        times (list[float]): List of time steps.
    """

    provide = SOP.TimeStepper

    @property
    def null_increment(self):
        return -1.0e150

    def __init__(self, times):
        super().__init__()
        self.current = 0
        self.level = 0
        assert sorted(times) == times
        self.times = list(times)

    def setInitialStep(self, time, prec=1e-6):
        """Update the sequence of times by keeping only the times greater
        than `time`.

        Arguments:
            time (float): All values smaller or equal to `time` are removed.
        """
        size = len(self.times)
        times = self.times
        index = next(i for i in range(size) if times[i] > time)
        self.times = times[index:]

    def setFinalStep(self, time, prec=1e-6):
        """Update the sequence of times by keeping only the times lower
        than `time`.

        Arguments:
            time (float): All values gretter or equal to `time` are removed.
        """
        size = len(self.times)
        times = self.times
        index = next(i for i in range(size) if times[i] >= time + prec)
        self.times = times[:index]

    def hasFinished(self):
        """Tell if there are steps to be computed.

        Returns:
            bool: *True* if there is no step to be computed, *False* otherwise.
        """
        return self.current > len(self.times) - 1

    def insertStep(self, index, time):
        """Inserts the step to given index.

        Arguments:
            index (int): index to insert nex step
            time (float): time value to insert.
        """

        self.times.insert(index, time)

    def getPrevious(self):
        """Returns the previous step calculated.

        Returns:
            float: Previous time value.
        """
        if self.hasFinished() or self.current - 1 < 0:
            raise IndexError

        return self.times[self.current - 1]

    def getNext(self):
        """Returns the next step to be calculated.

        Returns:
            float: Next time value.
        """
        if self.hasFinished():
            raise IndexError
        return self.times[self.current]

    def getIncrement(self):
        """Returns the increment to next step to be calculated.

        Returns:
            float: increment to the next time value.
        """
        if self.hasFinished():
            raise IndexError
        if self.current == 0:
            return self.null_increment

        return self.times[self.current] - self.times[self.current - 1]

    def completed(self):
        """Register the current step as completed successfully."""
        self.current += 1

    def split(self, nb_step):
        """Split last time step in nb_step uniform sub-steps

        Arguments:
            nb_step (int): Number of sub-steps.
        """
        # manage substepping...
        # for example, step back to the previous step
        assert nb_step > 1 and self.level < 20

        time_prev = self.times[self.current - 1]
        time_next = self.getNext()
        time_step = (time_next - time_prev) / float(nb_step)

        self.level += 1

        time_add = time_next
        for i in range(max(0, nb_step - 1)):
            time_add -= time_step
            self.insertStep(self.current, time_add)

    def raiseError(self, exc):
        """Raise an error executing the last step.

        Arguments:
            exc (Exception): Error to be raised.
        """
        # manage substepping...
        # for example, step back to the previous step
        self.current = max(self.current - 1, 0)
        raise exc


if __name__ == "__main__":
    stepper = TimeStepper([1.0, 2.0, 3.0])
    assert not stepper.hasFinished()
    step = stepper.getNext()
    assert step == 1.0, step
    step = stepper.getNext()
    assert step == 1.0, step
    stepper.completed()
    step = stepper.getNext()
    assert step == 2.0, step
    step = stepper.getNext()
    assert step == 2.0, step
    try:
        stepper.raiseError(ValueError())
        assert False
    except ValueError:
        pass
    assert not stepper.hasFinished()
    step = stepper.getNext()
    assert step == 1.0, step
    stepper.completed()
    step = stepper.getNext()
    assert step == 2.0, step
    stepper.completed()
    step = stepper.getNext()
    assert step == 3.0, step
    step = stepper.getNext()
    assert step == 3.0, step
