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

from ..Supervis import ConvergenceError
from ..Utilities import no_new_attributes, profile
from .solver_features import SolverFeature
from .solver_features import SolverOptions as SOP

import numpy as np


class LineSearch(SolverFeature):
    """Line search methods"""

    provide = SOP.LineSearch
    required_features = [SOP.ResidualComputation, SOP.PhysicalState]

    param = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, keywords):
        super().__init__()

        self.param = keywords
        print("PAR: ", self.param)

    def activated(self):
        """Return True if LineSearch is activated

        Returns:
            bool: True if LineSearch is activated.
        """

        return self.param is not None

    @profile
    @SolverFeature.check_once
    def solve(self, solution, scaling=1.0):
        """Apply linear search.

        Arguments:
            solution (FieldOnNodes): Displacement solution.

        Returns:
            FieldOnNodes: Accelerated solution by linear search.
        """

        if self.activated():
            method = self.param["METHODE"]
            if method not in ("CORDE"):
                raise NotImplementedError(method)

            def _f(rho, solution=solution, scaling=scaling):
                self.phys_state.primal_step += rho * solution
                # compute residual
                resiComp = self.get_feature(SOP.ResidualComputation)
                resi_state, varState, stressState = resiComp.computeResidual(scaling)
                self.phys_state.primal_step -= rho * solution
                return resi_state.resi.dot(solution), varState, stressState

            def _proj(rho, rhomin, rhomax, rhoexm, rhoexp):
                rhotmp = rho
                if rhotmp < rhomin:
                    rho = rhomin
                if rhotmp > rhomax:
                    rho = rhomax
                if rhotmp < 0.0 and rhotmp >= rhoexm:
                    rho = rhoexm
                if rhotmp >= 0.0 and rhotmp <= rhoexp:
                    rho = rhoexp

                return rho

            # retrieve args
            itemax = self.param["ITER_LINE_MAXI"]
            rhomin = self.param["RHO_MIN"]
            rhomax = self.param["RHO_MAX"]
            rtol = self.param["RESI_LINE_RELA"]
            rhoexm = self.param["RHO_EXCL"]
            rhoexp = self.param["RHO_EXCL"]

            # Implementing Secant Method
            rhom, rho = 0.0, 1.0
            rhoopt = rho
            fm, _, _ = _f(rhom, solution)  # minus residual is returned
            fopt = np.finfo("float64").max
            tiny = np.finfo("float64").tiny
            fcvg = abs(rtol * fm)

            for iter in range(itemax):
                try:
                    f, varState, stressState = _f(rho)
                except Exception as e:
                    # do we already have an rhoopt ?
                    if iter > 0:
                        return rhoopt * solution, varState, stressState
                    else:
                        raise e

                # keep best rho
                if abs(f) < fopt:
                    rhoopt = rho
                    fopt = abs(f)

                    # converged ?
                    if abs(f) < fcvg:
                        return rhoopt * solution, varState, stressState

                rhotmp = rho
                if abs(f - fm) > tiny:
                    rho = (f * rhom - fm * rho) / (f - fm)
                    rho = _proj(rho, rhomin, rhomax, rhoexm, rhoexp)
                elif f * (rho - rhom) * (f - fm) <= 0.0:
                    rho = rhomax
                else:
                    rho = rhomin

                print("Iteration-%d, rho2 = %0.6f and f(rho2) = %0.6f" % (iter, rho, f))
                rhom = rhotmp
                fm = f

            raise ConvergenceError()

        return solution, None, None
