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

from ...Utilities import no_new_attributes
from ..Basics import ContextMixin, KeywordsStore, ProblemTypeMixin
from ..Basics import ProblemType as PBT


class IterationsSolver(ABC, ContextMixin, ProblemTypeMixin):
    """Solves a step, loops on iterations."""

    class SubType(IntFlag):
        """Types of time integrators."""

        Unset = auto()
        Newton = auto()
        Snes = auto()
        Raspen = auto()

    solver_type = SubType.Unset

    logManager = None
    current_incr = current_matrix = None
    matr_update_incr = prediction = None
    __setattr__ = no_new_attributes(object.__setattr__)

    @classmethod
    def factory(cls, problem_type, keywords):
        """Factory that creates the appropriate object.

        Args:
            problem_type (ProblemType): Type of physical problem.
            keywords (dict): Part of user keywords.

        Returns:
            *BaseIntegrator*: A relevant *BaseIntegrator* object.
        """
        assert problem_type == PBT.MecaStat, f"unsupported type: {problem_type}"
        method = KeywordsStore(keywords).get("METHODE", default="NEWTON").capitalize()
        for kls in cls.__subclasses__():
            if kls.solver_type.name == method:
                return kls(keywords)

    def __init__(self):
        super().__init__()

    def initialize(self):
        """Initialize the object for the next step."""
        self.current_incr = 0
        self.current_matrix = None
        self.conv_manager.initialize()
        iter_glob = self.conv_manager.setdefault("ITER_GLOB_MAXI")
        iter_glob.minValue = 1

    def setParameters(self, param):
        """Assign parameters from user keywords.

        Arguments:
            param (dict) : user keywords.
        """
        self.param = param
        self.prediction = self._get("NEWTON", "PREDICTION") or self._get("NEWTON", "MATRICE")
        assert self.prediction in ("ELASTIQUE", "TANGENTE"), f"unsupported value: "
        self.matr_update_incr = self._get("NEWTON", "REAC_ITER", 1)

    def setLoggingManager(self, logManager):
        """Assign the logging manager.

        Arguments:
            logManager (LoggingManager): Logging manager.
        """
        self.logManager = logManager

    def update(self, primal_incr, resi_fields=None, callback=None):
        """Update the physical state.

        Arguments:
            primal_incr (FieldOnNodes): Displacement increment.
            resi_fields (dict of FieldOnNodes): Fields of residual values
        """

        self.phys_state.primal_step += primal_incr

        for key, field in resi_fields.items():
            self.phys_state.set(key, field)

        if callback:
            callback(primal_incr)

    def _resetMatrix(self, current_incr):
        """Reset matrix if needed

        Arguments:
            current_incr (int): index of the current increment

        """
        if (
            self.matr_update_incr > 0 and current_incr % self.matr_update_incr == 0
        ) or self.contact_manager:
            # make unavailable the current tangent matrix
            self.current_matrix = None

    def _setMatrixType(self, current_incr):
        """Set matrix type.

        Arguments:
            current_incr (int): index of the current increment

        Returns:
            str: Type of matrix to be computed.
        """
        if current_incr == 0:
            matrix_type = "PRED_" + self.prediction
        else:
            matrix_type = self._get("NEWTON", "MATRICE", "TANGENTE")

            self._resetMatrix(current_incr)
        return matrix_type

    @abstractmethod
    def solve(self, current_matrix, callback=None):
        """Solve a step.

        Raises:
            *ConvergenceError* exception in case of error.
        """
