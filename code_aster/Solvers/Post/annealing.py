# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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

from ...Objects import PostProcessing
from ..Basics import SolverOptions as SOP


class Annealing:
    """This object deals with annealing.

    This class can be used to calculate the mechanical state after annealing for
    non-linear mechanical calculations.

    It does:

    ..math::

        NewInternVar^{n} = f( InterVar^{n}, t^{n-1}, t^{n}, ExtVar^{n-1}, ExtVar^{n})
    """

    provide = SOP.PostStepHook
    required_features = [SOP.PhysicalProblem, SOP.PhysicalState]

    def __call__(self, nl_solver):
        previous = nl_solver.phys_state.getState(-1)
        current = nl_solver.phys_state
        post_process = PostProcessing(nl_solver.phys_pb)
        internVar_anneal = post_process.computeAnnealing(
            current.internVar,
            previous.time_curr,
            current.time_curr,
            previous.externVar,
            current.externVar,
        )
        current.set("VARI_ELGA", internVar_anneal)
        current.internVar = internVar_anneal
