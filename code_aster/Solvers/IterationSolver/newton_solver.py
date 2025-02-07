# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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

from ...Supervis import ConvergenceError
from ...Utilities import no_new_attributes, profile
from .convergence_manager import ConvergenceManager
from .incremental_solver import IncrementalSolver
from .iteration_solver import BaseInterationSolver
from .line_search import LineSearch


class NewtonSolver(BaseInterationSolver):
    """Solves a step, loops on iterations."""

    __needs__ = ("problem", "state", "keywords", "oper", "linear_solver", "contact")
    solver_type = BaseInterationSolver.SubType.Newton
    # FIXME: merge IncrementalSolver into NewtonSolver
    _incr_solv = _converg = _line_search = None
    __setattr__ = no_new_attributes(object.__setattr__)

    @classmethod
    def builder(cls, context):
        """Builder of RaspenSolver object.

        Args:
            context (Context): Context of the problem.

        Returns:
            instance: New object.
        """
        instance = cls()
        instance.context = context
        instance._post_init()
        instance._converg = ConvergenceManager.builder(context)
        instance._line_search = LineSearch.builder(context)
        instance._incr_solv = IncrementalSolver.builder(context)
        instance._incr_solv.share(instance._converg, instance._line_search)
        return instance

    def __init__(self):
        super().__init__()

    def initialize(self):
        """Initialize the object for the next step."""
        super().initialize()
        self._converg.initialize()
        iter_glob = self._converg.setdefault("ITER_GLOB_MAXI")
        iter_glob.minValue = 1

    def update(self, primal_incr, resi_fields=None, callback=None):
        """Update the physical state.

        Arguments:
            primal_incr (FieldOnNodes): Displacement increment.
            resi_fields (dict of FieldOnNodes): Fields of residual values
        """

        self.state.primal_step += primal_incr

        for key, field in resi_fields.items():
            self.state.set(key, field)

        if callback:
            callback(primal_incr)

    def _resetMatrix(self):
        """Reset matrix if needed

        Arguments:
            current_incr (int): index of the current increment
        """
        if (
            self._matr_update_incr > 0 and self.current_incr % self._matr_update_incr == 0
        ) or self.contact:
            # make unavailable the current tangent matrix
            self.current_matrix = None

    # @profile
    def solve(self, current_matrix, callback=None):
        """Solve a step.

        Raises:
            *ConvergenceError* exception in case of error.
        """
        self.current_matrix = current_matrix

        iter_glob = self._converg.setdefault("ITER_GLOB_MAXI")

        self.oper.initialize()
        while not self._converg.isFinished():
            iter_glob.value = self.current_incr

            # Select type of matrix
            matrix_type = self.matrix_type
            if self.current_incr > 0:
                self._resetMatrix()

            # Should the iteration be executed even if the solver converged ?
            force = self.oper.shouldExecuteIteration(self.current_incr)

            # Solve current iteration
            primal_incr, self.current_matrix, resi_fields = self._incr_solv.solve(
                matrix_type, self.current_matrix, force
            )

            # Update
            self.update(primal_incr, resi_fields, callback)

            if self.current_incr > 0:
                self.logManager.printConvTableRow(
                    [
                        self.current_incr - 1,
                        self._converg.get("RESI_GLOB_RELA"),
                        self._converg.get("RESI_GLOB_MAXI"),
                        self._converg.get("RESI_GEOM"),
                        matrix_type,
                    ]
                )
            self.current_incr += 1

        if not self._converg.isConverged():
            raise ConvergenceError("MECANONLINE9_7")

        self.oper.finalize()
        self._resetMatrix()
        return self.current_matrix
