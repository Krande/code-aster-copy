# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

# person_in_charge: mathieu.courtois@edf.fr
"""
:py:class:`LinearSolver` --- Base of linear solver objects
**********************************************************
"""

from ..Cata.Commons.c_solveur import C_SOLVEUR
from ..Objects import (
    GcpcSolver,
    LdltSolver,
    MultFrontSolver,
    MumpsSolver,
    PetscSolver,
    LinearSolver,
)
from ..Utilities import injector, logger


@injector(LinearSolver)
class ExtendedLinearSolver:
    @classmethod
    def factory(cls, mcf=None, **kwargs):
        """Create the solver object from the SOLVEUR factor keyword.

        Arguments:
            mcf (list|tuple|dict): Convenient option to pass `(kwargs,)`, the value
                returned for a factor keyword.
            kwargs (dict): Valid SOLVEUR keywords (syntax checked, with defaults).

        Returns:
            :class:`~code_aster.Objects.LinearSolver` (derivated of):
            Instance of a solver or *None* if none is selected.
        """
        if mcf:
            if isinstance(mcf, (list, tuple)):
                mcf = mcf[0]
            if isinstance(mcf, dict):
                kwargs = mcf
        name = kwargs.get("METHODE")
        klass = None
        for sub in cls.__subclasses__():
            if sub.solverName() == name:
                klass = sub
                break
        assert klass, f"Unknown solver: {name}"
        solver = klass(**kwargs)
        return solver


class LinearSolverExt:
    """Base object for LinearSolver."""

    cata_sdj = "SD.sd_solveur.sd_solveur"
    _name = _init = None

    def __init__(self, *args, **kwargs):
        assert len(args) <= 1, "at most one argument is expected"
        self._init(*args)
        mcf = C_SOLVEUR("MECA_STATIQUE")
        keywords = dict(METHODE=self._name)
        mcf.addDefaultKeywords(keywords)
        keywords.update(kwargs)
        logger.debug("solver init with:", keywords)
        self.setKeywords(keywords)

    @classmethod
    def solverName(cls):
        """Return the solver name.

        Returns:
            str: Solver name, matching METHODE value.
        """
        return cls._name

@injector(GcpcSolver)
class ExtendedGcpcSolver(LinearSolverExt):
    _name = "GCPC"
    _init = GcpcSolver.__init__


@injector(LdltSolver)
class ExtendedLdltSolver(LinearSolverExt):
    _name = "LDLT"
    _init = LdltSolver.__init__


@injector(MultFrontSolver)
class ExtendedMultFrontSolver(LinearSolverExt):
    _name = "MULT_FRONT"
    _init = MultFrontSolver.__init__


@injector(MumpsSolver)
class ExtendedMumpsSolver(LinearSolverExt):
    _name = "MUMPS"
    _init = MumpsSolver.__init__


@injector(PetscSolver)
class ExtendedPetscSolver(LinearSolverExt):
    _name = "PETSC"
    _init = PetscSolver.__init__
