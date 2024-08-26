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

from ...Supervis import ConvergenceError
from ...Utilities import no_new_attributes, profile, logger
from ..Basics import SolverFeature
from ..Basics import SolverOptions as SOP

import numpy as np


class LineSearch(SolverFeature):
    """Line search methods"""

    provide = SOP.LineSearch

    required_features = [SOP.PhysicalProblem, SOP.PhysicalState]

    optional_features = [SOP.OperatorsManager]

    param = solve = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, keywords):
        super().__init__()

        self.param = keywords
        self.solve = None

    def activated(self):
        """Return True if LineSearch is activated

        Returns:
            bool: True if LineSearch is activated.
        """

        return self.param is not None and self.param["ITER_LINE_MAXI"] > 0

    def setup(self):
        """Set up the line search object"""
        self.solve = self.__solve

    @profile
    @SolverFeature.check_once
    def __solve(self, solution, scaling=1.0):
        """Apply linear search for mechanical problems.

        Arguments:
            solution (FieldOnNodes): Displacement solution.

        Returns:
            FieldOnNodes: Accelerated solution by linear search.
        """

        if self.activated():
            method = self.param["METHODE"]

            def _f(rho, solution=solution, scaling=scaling):
                self.phys_state.primal_step += rho * solution
                # compute residual
                opersManager = self.get_feature(SOP.OperatorsManager)
                resi_state = opersManager.getResidual(scaling)
                self.phys_state.primal_step -= rho * solution
                return -resi_state.resi.dot(solution)

            def _proj(rho, param=self.param):
                rhotmp = rho
                if rhotmp < param["RHO_MIN"]:
                    rho = param["RHO_MIN"]
                if rhotmp > param["RHO_MAX"]:
                    rho = param["RHO_MAX"]
                if rhotmp < 0.0 and rhotmp >= -param["RHO_EXCL"]:
                    rho = -param["RHO_EXCL"]
                if rhotmp >= 0.0 and rhotmp <= param["RHO_EXCL"]:
                    rho = param["RHO_EXCL"]

                return rho

            # retrieve args
            f0 = _f(0.0, solution)
            fopt = np.finfo("float64").max
            tiny = np.finfo("float64").tiny
            fcvg = abs(self.param["RESI_LINE_RELA"] * f0)
            iteropt = -1

            assert method in ("CORDE", "MIXTE"), method
            rhom, rho = 0.0, 1.0
            rhoopt = rho
            fm = f0
            if method == "MIXTE":
                if f0 <= 0.0:
                    sens = 1.0
                else:
                    sens = -1.0
                rhoneg = 0.0
                fneg = sens * f0
                bpos = False
            else:
                sens = 1.0

            for iter in range(self.param["ITER_LINE_MAXI"] + 1):
                try:
                    f = _f(sens * rho)
                except Exception as e:
                    # do we already have an rhoopt ?
                    if iter > 0:
                        return rhoopt * solution
                    else:
                        raise e
                # keep best rho
                # zbopti
                if abs(f) <= fopt:
                    rhoopt = rho
                    fopt = abs(f)
                    iteropt = iter
                    # converged ?
                    if abs(f) < fcvg:
                        logger.info(
                            "Linesearch: iter = %d, rho = %0.6f and f(rho) = %0.6f" % (iter, rho, f)
                        )
                        return rhoopt * solution

                rhotmp = rho
                if method == "CORDE":
                    if abs(f - fm) > tiny:
                        rho = (f * rhom - fm * rho) / (f - fm)
                        rho = _proj(rho)
                    elif f * (rho - rhom) * (f - fm) <= 0.0:
                        rho = self.param["RHO_MAX"]
                    else:
                        rho = self.param["RHO_MIN"]
                elif method == "MIXTE":
                    # zbborn
                    if np.sign(f) == np.sign(f0):
                        rhoneg = rho
                    else:
                        rhopos = rho
                        bpos = True
                    if not bpos:
                        rhom = rho
                        rho = 3 * rhom
                    else:
                        # zbroot
                        if abs(f) >= abs(fm):
                            # en cas de non pertinence des iteres : dichotomie
                            rho = 0.5 * (rhoneg + rhopos)
                        else:
                            # interpolation lineaire
                            if abs(rho - rhom) > tiny:
                                p1 = (f - fm) / (rho - rhom)
                                p0 = fm - p1 * rhom
                                if abs(p1) <= abs(fm) / (rhopos + rhom):
                                    rho = 0.5 * (rhoneg + rhopos)
                                else:
                                    rho = -p0 / p1
                            else:
                                logger.info(
                                    "Linesearch: iter = %d, rho = %0.6f and f(rho) = %0.6f"
                                    % (iteropt, rhoopt, fopt)
                                )
                                return rhoopt * solution
                    # zbproj
                    if rho < rhoneg:
                        if bpos:
                            print("ici ici")
                            rho = 0.5 * (rhoneg + rhopos)
                        else:
                            logger.info(
                                "LinesearchB: iter = %d, rho = %0.6f and f(rho) = %0.6f"
                                % (iteropt, rhoopt, fopt)
                            )
                            return rhoopt * solution
                    if bpos and rho > rhopos:
                        rho = 0.5 * (rhoneg + rhopos)
                    # zbinte
                    rho = sens * _proj(sens * rho)
                rhom = rhotmp
                fm = f
            logger.info(
                "Linesearch: iter = %d, rho = %0.6f and f(rho) = %0.6f" % (iteropt, rhoopt, fopt)
            )
            return rhoopt * solution

        return solution
