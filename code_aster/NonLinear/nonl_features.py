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
Objects used to build non linear operators.
"""

from .base_features import BaseFeature, BaseFeaturesOptions


class NonLinearOptions(BaseFeaturesOptions):
    """Enumeration of non linear options."""

    # description of the system (or study, problem)
    PhysicalProblem = 0x0001
    PhysicalState = 0x0002
    # object that store the results
    Storage = 0x0004

    ProblemSolver = 0x0008
    TimeStepper = 0x0010
    StepSolver = 0x0020
    ConvergenceCriteria = 0x040
    ConvergenceManager = 0x080
    LinearSolver = 0x0100
    IncrementalSolver = 0x0200
    Contact = 0x0400


class NonLinearFeature(BaseFeature):
    """Feature object for non linear operators."""

    options = NonLinearOptions

    # convenient shortcuts properties
    @property
    def phys_pb(self):
        """PhysicalProblem: current problem description."""
        return self.get_feature(NonLinearOptions.PhysicalProblem)

    @property
    def phys_state(self):
        """PhysicalState: current state."""
        return self.get_feature(NonLinearOptions.PhysicalState)
