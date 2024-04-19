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
from ...Utilities import no_new_attributes, profile
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
        if self.phys_pb.isMechanical():
            self.solve = self.__solve_meca
        elif self.phys_pb.isThermal():
            self.solve = self.__solve_ther

    @profile
    @SolverFeature.check_once
    def __solve_meca(self, solution, scaling=1.0):
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
                return resi_state.resi.dot(solution)

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
            if method == "CORDE":
                rhom, rho = 0.0, 1.0
                rhoopt = rho
                fm = f0
            elif method == "MIXTE":
                raise NotImplementedError()

            for iter in range(self.param["ITER_LINE_MAXI"] + 1):
                try:
                    f = _f(rho)
                except Exception as e:
                    # do we already have an rhoopt ?
                    if iter > 0:
                        return rhoopt * solution
                    else:
                        raise e
                # keep best rho
                if abs(f) <= fopt:
                    rhoopt = rho
                    fopt = abs(f)
                    iteropt = iter
                    # converged ?
                    if abs(f) < fcvg:
                        print(
                            "Linesearch: iter = %d, rho = %0.6f and f(rho) = %0.6f" % (iter, rho, f)
                        )
                        return rhoopt * solution

                rhotmp = rho
                if abs(f - fm) > tiny:
                    rho = (f * rhom - fm * rho) / (f - fm)
                    rho = _proj(rho)
                elif f * (rho - rhom) * (f - fm) <= 0.0:
                    rho = self.param["RHO_MAX"]
                else:
                    rho = self.param["RHO_MAX"]

                # print("Iteration-%d, rho2 = %0.6f and f(rho2) = %0.6f" % (iter, rho, f))
                rhom = rhotmp
                fm = f
            print("Linesearch: iter = %d, rho = %0.6f and f(rho) = %0.6f" % (iteropt, rhoopt, fopt))
            return rhoopt * solution

        return solution

    @profile
    @SolverFeature.check_once
    def __solve_ther(self, solution, scaling=1.0):
        """Apply linear search for thermal problems.

        Arguments:
            solution (FieldOnNodes): Temperature solution.

        Returns:
            FieldOnNodes: Accelerated solution by linear search.
        """

        if self.activated():

            def _f(rho, solution=solution, scaling=scaling):
                self.phys_state.primal_step += rho * solution
                opersManager = self.get_feature(SOP.OperatorsManager)
                resi_state = opersManager.getResidual(scaling)
                self.phys_state.primal_step -= rho * solution
                testm = resi_state.resi.norm("NORM_INFINITY")
                return resi_state.resi.dot(solution), testm

            # retrieve args
            f0, testm = _f(0.0, solution)

            rho0 = 0.0
            rho = 1.0

            for iter in range(self.param["ITER_LINE_MAXI"] + 1):

                f1, testm = _f(rho)

                if testm < self.param["RESI_LINE_RELA"]:
                    return rho * solution

                if iter == 0:
                    ffinal = f1
                    rhof = 1.0

                if abs(f1) < abs(ffinal):
                    ffinal = f1
                    rhof = rho

                rhot = rho

                if abs(f1 - f0) > np.finfo("float64").tiny:
                    rho = -(f0 * rhot - f1 * rho0) / (f1 - f0)
                    # print(
                    #     "Linesearch: f1 = {}, rho0 = {}, f0 = {}, rhot = {}, rho = {}".format(f1, rho0, f0, rhot, rho),
                    #     flush=True
                    # )
                    if rho < self.param["RHO_MIN"]:
                        rho = self.param["RHO_MIN"]
                    if rho > self.param["RHO_MAX"]:
                        rho = self.param["RHO_MAX"]
                    if abs(rho - rhot) < 1.0e-8:
                        break
                else:
                    break

                rho0 = rhot
                f0 = f1

            rho = rhof
            f1 = ffinal

            return rho * solution

        return solution
