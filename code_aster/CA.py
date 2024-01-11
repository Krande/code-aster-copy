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

# person_in_charge: mathieu.courtois@edf.fr

"""
This modules provides a compatibility module to emulate
"""

from .Commands import *
from .Commands.debut import init
from .Commands.fin import FIN as close
from .Objects import *
from .ObjectsExt import DataStructure
from .Supervis import (
    AsterError,
    ContactError,
    ConvergenceError,
    IntegrationError,
    SolverError,
    TimeLimitError,
    saveObjects,
)
from .Utilities import MPI, TestCase
from .Utilities.version import get_version
from .Utilities.version import version_info

# may happen when building the doc
__version__ = get_version() if version_info else ""

del get_version, version_info
