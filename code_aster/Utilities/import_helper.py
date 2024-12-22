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

"""
This module gives wrapper objects to use PETSc and medcoupling as optional modules.
The modules are actually imported at the first use.

Please *always* import PETSc or medcoupling using this wrapper rather than directly
by doing: ``from ..Utilities import PETSc`` or ``from ..Utilities import medcoupling``.
"""

# aslint: disable=C4008


class _PETScMeta(type):
    """Meta class for petsc4py wrapping."""

    _init = False
    _mod = None

    def __getattr__(cls, attr):
        if not cls._init:
            import petsc4py
            from petsc4py import PETSc as origin

            cls._init = True
            cls._mod = origin
        return getattr(cls._mod, attr)


class _SLEPcMeta(type):
    """Meta class for slepc4py wrapping."""

    _init = False
    _mod = None

    def __getattr__(cls, attr):
        if not cls._init:
            import slepc4py
            from slepc4py import SLEPc as origin

            cls._init = True
            cls._mod = origin
        return getattr(cls._mod, attr)


class PETSc(metaclass=_PETScMeta):
    """Wrapper to petsc4py.PETSc"""


class SLEPc(metaclass=_SLEPcMeta):
    """Wrapper to slepc4py.SLEPc"""


class _medcouplingMeta(type):
    """Meta class for medcoupling module wrapping."""

    _init = False
    _mod = None

    def __getattr__(cls, attr):
        if not cls._init:
            import medcoupling as origin

            cls._init = True
            cls._mod = origin
        return getattr(cls._mod, attr)


class _ParaMEDMEMMeta(type):
    """Meta class for ParaMEDMEM module wrapping."""

    _init = False
    _mod = None

    def __getattr__(cls, attr):
        if not cls._init:
            import ParaMEDMEM as origin

            cls._init = True
            cls._mod = origin
        return getattr(cls._mod, attr)


class medcoupling(metaclass=_medcouplingMeta):
    """Wrapper to medcoupling"""


class ParaMEDMEM(metaclass=_ParaMEDMEMMeta):
    """Wrapper to ParaMEDMEM"""
