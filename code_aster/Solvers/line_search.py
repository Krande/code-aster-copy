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

from ..Utilities import no_new_attributes, profile
from .solver_features import SolverFeature
from .solver_features import SolverOptions as SOP


class LineSearch(SolverFeature):
    """Line search methods"""

    provide = SOP.LineSearch
    param = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, keywords):
        super().__init__()

        self.param = keywords
        print("PAR: ", self.param)

    @profile
    @SolverFeature.check_once
    def solve(self, field):
        """Apply linear search.

        Arguments:
            field (FieldOnNodes): Displacement field.

        Returns:
            FieldOnNodes: Accelerated field by linear search.
        """

        if self.param:

            def _f(rho, field=field, scaling=scaling):
                self.phys_state.primal_step += rho * field
                resi_state, _, _ = self.computeResidual(scaling)
                self.phys_state.primal_step -= rho * field

            return -resi_state.resi.dot(field)

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

            # retrieve args
            itemax = self.param["ITER_LINE_MAXI"]
            rhomin = self.param["RHO_MIN"]
            rhomax = self.param["RHO_MAX"]
            rtol = self.param["RESI_LINE_RELA"]
            rhoexm = self.param["RHO_EXCL"]
            rhoexp = self.param["RHO_EXCL"]

            # Implementing Secant Method
            rhom, rho = 0.0, 1.0
            step = 1
            fm = -residual.dot(field)  # minus residual is returned
            itemax = 30
            fopt = np.finfo("float64").max
            tiny = np.finfo("float64").tiny
            rhoopt = 1.0
            fcvg = abs(rtol * fm)

            while True:
                try:
                    f = _f(rho)
                except Exception as e:
                    # do we already have an rhoopt ?
                    if step > 1:
                        return rhoopt * field
                    else:
                        raise e
                # keep best rho
                if abs(f) < fopt:
                    rhoopt = rho
                    fopt = abs(f)
                    # converged ?
                    if abs(f) < fcvg:
                        return rhoopt * field
                rhotmp = rho
                if abs(f - fm) > tiny:
                    rho = (f * rhom - fm * rho) / (f - fm)
                    _proj(rho, rhomin, rhomax, rhoexm, rhoexp)
                elif f * (rho - rhom) * (f - fm) <= 0.0:
                    rho = rhomax
                else:
                    rho = rhomin

                print("Iteration-%d, rho2 = %0.6f and f(rho2) = %0.6f" % (step, rho, f))
                rhom = rhotmp
                fm = f
                step = step + 1

                if step > itemax:
                    raise RuntimeError("No convergence in line search")

        return field
