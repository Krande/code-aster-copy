# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

# person_in_charge: nicolas.pignet@edf.fr

"""
This module gives common utilities for MPI.

Need only mpi4py package
"""

from collections import defaultdict
from .base_utils import ReadOnlyDict

try:
    from .aster_config import config as _cfg
    config = ReadOnlyDict(**_cfg)
    del _cfg
except ImportError:
    config = defaultdict(lambda: None)

def haveMPI():
    """Tell if the library is built with MPI support.

    Returns:
    bool: *True* if use MPI librairies, *False* else
    """
    return config.get('ASTER_HAVE_MPI', 0) == 1

try:
    from mpi4py import MPI
except:
    class FAKE_COMM_WORLD:
        """
        This class FAKE_COMM_WORLD contains methods for compatibility with sequential libraries

        Use MPI.COMM_WORLD as mpi4py to use mpi methods

        Some methods can be missing (add them here with the same name and arguments than mpi4py)
        """

        def __init__(self):
            if haveMPI():
                raise RuntimeError("mpi4py is mandatory for mpi execution")

        def Get_rank(self):
            return 0

        def Get_size(self):
            return 1

        def Barrier(self):
            return

        def bcast(self, data, root):
            return data

        def Bcast(self, data, root):
            return data

        def gather(self, data, root):
            return data

    class MPI:
        """
        This class MPI is an encapsulation of mpi4py for sequential libraries.

        The same API than mpi4py is used

        It is equivalent to 'from mpi4py import MPI' of parallel libraries
        """

        COMM_WORLD = FAKE_COMM_WORLD()
