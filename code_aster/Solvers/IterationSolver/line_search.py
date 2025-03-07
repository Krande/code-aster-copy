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

from abc import ABC, abstractmethod
from enum import IntFlag, auto

import numpy as np

from ...Utilities import logger, no_new_attributes, profile
from ..Basics import ContextMixin


class LineSearchType(IntFlag):
    """Line Search types."""

    Unset = auto()
    Corde = auto()
    Mixte = auto()
    Pilotage = auto()

    @classmethod
    def by_name(cls, name):
        """Return an option value by its name.
        *AttributeError* is raised if the option does not exist.

        Arguments:
            name (str): Option name.

        Returns:
            int: Option value.
        """
        return getattr(cls, name)


class BaseLineSearch(ABC, ContextMixin):
    """Base for line search methods"""

    __needs__ = ("keywords", "state", "oper")
    linesearch_type = LineSearchType.Unset
    _cached_limits = None
    __setattr__ = no_new_attributes(object.__setattr__)

    @classmethod
    def factory(cls, context):
        """Factory that creates the appropriate object.

        Args:
            context (Context): Context of the problem.

        Returns:
            *BaseLineSearch*: New object.
        """
        method = context.get_keyword("RECH_LINEAIRE", "METHODE", "CORDE").capitalize()
        searched = LineSearchType.by_name(method)
        for kls in cls.__subclasses__():
            if searched in kls.linesearch_type:
                return kls.builder(context)
        raise TypeError(f"no candidate for cls={cls}, method: {method}")

    def __init__(self):
        super().__init__()
        self._cached_limits = None

    def _get(self, keyword=None, default=None):
        return self.get_keyword("RECH_LINEAIRE", keyword, default)

    def isEnabled(self):
        """bool: *True* if LineSearch is activated."""
        return self._get("ITER_LINE_MAXI", 0) > 0

    def compute_f(self, rho, solution, scaling=1.0):
        """Compute functional"""
        self.state.primal_step += rho * solution
        # compute residual
        resi_state = self.oper.getResidual(scaling)
        self.state.primal_step -= rho * solution
        return -resi_state.resi.dot(solution)

    def check_limits(self, rho):
        if self._cached_limits is None:
            self._cached_limits = [self._get(i) for i in ("RHO_MIN", "RHO_MAX", "RHO_EXCL")]
        min, max, excl = self._cached_limits
        if rho < min:
            rho = min
        if rho > max:
            rho = max
        if -excl <= rho < 0.0:
            rho = -excl
        if 0.0 <= rho <= excl:
            rho = excl
        return rho

    @abstractmethod
    def solve(self, solution, scaling=1.0):
        """Apply linear search for mechanical problems.

        Arguments:
            solution (FieldOnNodes): Displacement solution.
            scaling (float): Scaling factor for Lagrange multipliers (default: 1.0).

        Returns:
            FieldOnNodes: Accelerated solution by linear search.
        """


class SecantLineSearch(BaseLineSearch):
    """Line search for CORDE & MIXTE methods."""

    linesearch_type = LineSearchType.Corde | LineSearchType.Mixte
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        super().__init__()

    @profile
    def solve(self, solution, scaling=1.0):
        """Apply linear search for mechanical problems.

        Arguments:
            solution (FieldOnNodes): Displacement solution.
            scaling (float): Scaling factor for Lagrange multipliers (default: 1.0).

        Returns:
            FieldOnNodes: Accelerated solution by linear search.
        """
        if not self.isEnabled():
            return solution

        method = self._get("METHODE")
        assert method in ("CORDE", "MIXTE"), method
        f0 = self.compute_f(0.0, solution)
        fopt = np.finfo("float64").max
        tiny = np.finfo("float64").tiny
        fcvg = abs(self._get("RESI_LINE_RELA") * f0)
        iteropt = -1

        rhom, rho = 0.0, 1.0
        rhoopt = rho
        fm = f0
        if method == "MIXTE":
            if f0 <= 0.0:
                sens = 1.0
            else:
                sens = -1.0
            rhoneg = 0.0
            bpos = False
        else:
            sens = 1.0

        for iter in range(self._get("ITER_LINE_MAXI") + 1):
            try:
                f = self.compute_f(sens * rho, solution)
            except Exception:
                # do we already have an rhoopt ?
                if iter > 0:
                    return rhoopt * solution
                raise
            # keep best rho
            # zbopti
            if abs(f) <= fopt:
                rhoopt = rho
                fopt = abs(f)
                iteropt = iter
                # converged ?
                if abs(f) < fcvg:
                    logger.debug("Linesearch: iter = %d, rho = %0.6f, f(rho) = %0.6f", iter, rho, f)
                    return rhoopt * solution

            rhotmp = rho
            if method == "CORDE":
                if abs(f - fm) > tiny:
                    rho = (f * rhom - fm * rho) / (f - fm)
                    rho = self.check_limits(rho)
                elif f * (rho - rhom) * (f - fm) <= 0.0:
                    rho = self._get("RHO_MAX")
                else:
                    rho = self._get("RHO_MIN")
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
                            logger.debug(
                                "Linesearch: iter = %d, rho = %0.6f, f(rho) = %0.6f",
                                iteropt,
                                rhoopt,
                                fopt,
                            )
                            return rhoopt * solution
                # zbproj
                if rho < rhoneg:
                    if bpos:
                        rho = 0.5 * (rhoneg + rhopos)
                    else:
                        logger.debug(
                            "Linesearch: iter = %d, rho = %0.6f, f(rho) = %0.6f",
                            iteropt,
                            rhoopt,
                            fopt,
                        )
                        return rhoopt * solution
                if bpos and rho > rhopos:
                    rho = 0.5 * (rhoneg + rhopos)
                # zbinte
                rho = sens * self.check_limits(sens * rho)
            rhom = rhotmp
            fm = f
        logger.debug("Linesearch: iter = %d, rho = %0.6f, f(rho) = %0.6f", iteropt, rhoopt, fopt)
        return rhoopt * solution
