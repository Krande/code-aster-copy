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

from ...Objects import HHO, PostProcessing
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

    def __init__(self) -> None:
        self._enabled = None

    def __call__(self, nl_solver):
        if self._enabled is None:
            self._enabled = nl_solver.phys_pb.getBehaviourProperty().hasAnnealing()
        if not self._enabled:
            return

        try:
            previous = nl_solver.phys_state.getState(-1)
        except IndexError:
            # No post-pro at initial step
            return
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


class ComputeHydr:
    """Hook to compute HYDR_ELGA."""

    provide = SOP.PostStepHook

    def __init__(self) -> None:
        self._enabled = None

    def __call__(self, nl_solver):
        if self._enabled is None:
            self._enabled = nl_solver.phys_pb.getBehaviourProperty().hasBehaviour("THER_HYDR")
        if not self._enabled:
            return

        current = nl_solver.phys_state
        post = PostProcessing(nl_solver.phys_pb)
        try:
            hydr_prev = current.getState(-1).auxiliary["HYDR_ELGA"]
            hydr_curr = post.computeHydration(
                current.primal_prev,
                current.primal_curr,
                current.time_prev,
                current.time_curr,
                hydr_prev,
            )
        except IndexError:
            hydr_curr = current.createFieldOnCells(nl_solver.phys_pb, "ELGA", "HYDR_R")
        current.set("HYDR_ELGA", hydr_curr)


class PostHHO:
    """Compute the true primal field from HHO unknowns."""

    provide = SOP.PostStepHook
    _field_name = None

    def __init__(self) -> None:
        self._enabled = None
        self._hho = None

    def __call__(self, nl_solver):
        """Hook to compute HHO_DEPL"""
        if self._enabled is None:
            self._enabled = nl_solver.phys_pb.getModel().existsHHO()
            self._hho = HHO(nl_solver.phys_pb)
        if not self._enabled:
            return

        current = nl_solver.phys_state
        hho_field = self._hho.projectOnLagrangeSpace(current.primal_curr)
        current.set(self._field_name, hho_field)


class ComputeDisplFromHHO(PostHHO):
    """Compute the displacement field from HHO unknowns."""

    _field_name = "HHO_DEPL"


class ComputeTempFromHHO(PostHHO):
    """Compute the temperature field from HHO unknowns."""

    _field_name = "HHO_TEMP"
