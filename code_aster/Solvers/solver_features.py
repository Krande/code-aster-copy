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

"""
Objects used to build generic problems solvers.
"""

from enum import IntFlag, auto
from .base_features import BaseFeature


class SolverOptions(IntFlag):
    """Enumeration of non linear options."""

    Nothing = 0
    # main object
    ProblemSolver = auto()
    # description of the system (or study, problem)
    PhysicalProblem = auto()
    PhysicalState = auto()
    # storage of the results
    Storage = auto()
    # time steps management
    TimeStepper = auto()
    # solves a step
    StepSolver = auto()
    # criteria that must be checked for a step
    ConvergenceCriteria = auto()
    ConvergenceManager = auto()
    # solves an increment
    IncrementalSolver = auto()
    # linear solver
    LinearSolver = auto()
    # line serach
    LineSearch = auto()
    # contact object
    Contact = auto()
    # container of keywords
    Keywords = auto()

    # flag for notification support
    EventSource = auto()

    # hooks called after a step
    PostStepHook = auto()


class SolverFeature(BaseFeature):
    """Feature object for non linear operators."""

    # convenient shortcuts properties
    @property
    def phys_pb(self):
        """PhysicalProblem: current problem description."""
        return self.get_feature(SolverOptions.PhysicalProblem)

    @property
    def phys_state(self):
        """PhysicalState: current state."""
        return self.get_feature(SolverOptions.PhysicalState)
