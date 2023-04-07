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

from .base_features import BaseFeature, BaseFeaturesOptions


class SolverOptions(BaseFeaturesOptions):
    """Enumeration of non linear options."""

    # main object
    ProblemSolver = 0x1
    # description of the system (or study, problem)
    PhysicalProblem = 0x2
    PhysicalState = 0x4
    # storage of the results
    Storage = 0x8
    # time steps management
    TimeStepper = 0x10
    # solves a step
    StepSolver = 0x20
    # criteria that must be checked for a step
    ConvergenceCriteria = 0x40
    ConvergenceManager = 0x80
    # solves an increment
    IncrementalSolver = 0x100
    # linear solver
    LinearSolver = 0x200
    # contact object
    Contact = 0x400
    # container of keywords
    Keywords = 0x800

    # flag added "for a step" object
    ForStep = 0x1000
    # flag added "for an increment" object
    ForIncr = 0x2000


class SolverFeature(BaseFeature):
    """Feature object for non linear operators."""

    options = SolverOptions

    # convenient shortcuts properties
    @property
    def phys_pb(self):
        """PhysicalProblem: current problem description."""
        return self.get_feature(SolverOptions.PhysicalProblem)

    @property
    def phys_state(self):
        """PhysicalState: current state."""
        return self.get_feature(SolverOptions.PhysicalState)
