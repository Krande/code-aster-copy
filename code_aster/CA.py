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
This module provides a low level interface to the most of the code_aster objects.

.. code-block:: python

    >>> import code_aster
    >>> # 'code_aster.rc' object is available for customization.
    >>> # the following line automatically starts 'init()':
    >>> from code_aster import CA
    >>> # 'CA' provides all user objects and commands.

In a standard commands file, without advanced Python usage:

.. code-block:: python

    >>> import code_aster
    >>> from code_aster.Commands import *


"""

# from .version import get_version, get_version_desc

import atexit
from .Utilities.rc import rc

if rc.initialize is None:
    rc.initialize = True
if rc.finalize is None:
    rc.finalize = rc.initialize

from .Commands import *
from .Commands.debut import init
from .Commands.fin import close
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
from .Utilities.version import __version__

if rc.initialize:
    init()

if rc.finalize:
    atexit.register(close)
