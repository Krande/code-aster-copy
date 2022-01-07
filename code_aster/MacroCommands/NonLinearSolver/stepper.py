# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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


class TimeStepper:
    """"Basic time stepper.

    Arguments:
        times (list[float]): List of time steps.
    """

    def __init__(self, times):
        self.current = 0
        self.times = times

    def updateTimes(self, time):
        """Update the sequence of times by keeping only the times greater
        than `time`.

        Arguments:
            time (float): All values smaller or equal to `time` are removed.
        """
        size = len(self.times)
        times = self.times
        index = next(i for i in range(size) if times[i] > time)
        self.times = times[index:]

    def hasFinished(self):
        """Tell if there are steps to be computed.

        Returns:
            bool: *True* if there is no step to be computed, *False* otherwise.
        """
        return self.current > len(self.times) - 1

    def getNext(self):
        """Returns the next step to be calculated.

        Returns:
            float: Next time value.
        """
        if self.hasFinished():
            raise IndexError
        return self.times[self.current]

    def completed(self):
        """Register the current step as completed successfully."""
        self.current += 1

    def raiseError(self, exc):
        """Raise an error executing the last step.

        Arguments:
            exc (Exception): Error to be raised.
        """
        # manage substepping...
        # for example, step back to the previous step
        self.current = max (self.current - 1, 0)
        raise exc


if __name__ == "__main__":
    stepper = TimeStepper([1., 2., 3.])
    assert not stepper.hasFinished()
    step = stepper.getNext()
    assert step == 1., step
    step = stepper.getNext()
    assert step == 1., step
    stepper.completed()
    step = stepper.getNext()
    assert step == 2., step
    step = stepper.getNext()
    assert step == 2., step
    try:
        stepper.raiseError(ValueError())
        assert False
    except ValueError:
        pass
    assert not stepper.hasFinished()
    step = stepper.getNext()
    assert step == 1., step
    stepper.completed()
    step = stepper.getNext()
    assert step == 2., step
    stepper.completed()
    step = stepper.getNext()
    assert step == 3., step
    step = stepper.getNext()
    assert step == 3., step
