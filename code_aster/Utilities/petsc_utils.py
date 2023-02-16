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

# person_in_charge: nicolas.tardieu@edf.fr

"""
This module gives common utilities for PETSc.
Please *always* import PETSc using this wrapper rather than directly.
"""

from libaster import _petscInitializeWithOptions, petscFinalize


class _PETScMeta(type):
    """Meta class for petsc4py wrapping."""

    _init = False
    _mod = None

    def __getattr__(cls, attr):
        if not cls._init:
            import petsc4py
            from petsc4py import PETSc

            cls._init = True
            cls._mod = PETSc
        return getattr(cls._mod, attr)


class PETSc(metaclass=_PETScMeta):
    """Wrapper to petsc4py.PETSc"""


def petscInitialize(options=""):
    """Starts the PETSc interface with options.

    Arguments:
        options[str]: PETSc options
    """
    _petscInitializeWithOptions(options)
