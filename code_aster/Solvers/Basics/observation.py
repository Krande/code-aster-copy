# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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

from ...Utilities import logger, no_new_attributes
from .bases import EventId, Observer
from .context import ContextMixin


class Observation(ContextMixin, Observer):
    """This object extracts values during a nonlinear process.

    - It contains the definition of the observable values:
      - what (field, components)
      - where (group of cells, nodes)
      - how (extraction, minimal, maximal, mean value)
      - when (during iterations for printing, to valid a time step, to be
        added in a table)
    """

    class Observable:
        field = components = location = operation = None
        target = None

    __needs__ = ("problem", "state")
    _data = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        self._data = []

    def notify(self, event):
        """Receive notification from event, store the data from a future use
        by actions.

        Arguments:
            event (EventSource): Object that sends the notification.
        """
        if event.eid & EventId.IterationSolver:
            logger.debug("+ check for observations during iterations")
            if self._data:  # loop on observables...
                self.state.debugPrint(recursive=True)
        elif event.eid & EventId.TimeStepper:
            logger.debug("+ check for observations for time stepper")
        elif event.eid & EventId.NonLinearOperator:
            logger.debug("+ check for observations at convergence")
        else:
            raise TypeError(f"unsupported event: eid={event.eid!r}")
