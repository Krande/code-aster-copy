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

# person_in_charge: mathieu.courtois@edf.fr
"""
Objects defined depending on the configuration
**********************************************

The availability of the objects depends on the configuration.
Some of them are not available with a sequential version.
Some others are only defined if a prerequisite is well configured.
"""

from .datastructure_py import UnavailableObject


try:
    from libaster import Mesh, ParallelMesh
except ImportError:

    class ParallelMesh(UnavailableObject):
        pass


try:
    from libaster import ConnectionMesh
except ImportError:

    class ConnectionMesh(UnavailableObject):
        pass


try:
    from libaster import CommGraph
except ImportError:

    class CommGraph(UnavailableObject):
        pass


try:
    from libaster import IncompleteMesh
except ImportError:

    class IncompleteMesh(UnavailableObject):
        pass


try:
    from libaster import ParallelFiniteElementDescriptor
except ImportError:

    class ParallelFiniteElementDescriptor(UnavailableObject):
        pass


try:
    from libaster import ParallelMechanicalLoadReal

except ImportError:

    class ParallelMechanicalLoadReal(UnavailableObject):
        pass


try:
    from libaster import ParallelMechanicalLoadFunction

except ImportError:

    class ParallelMechanicalLoadFunction(UnavailableObject):
        pass


try:
    from libaster import ParallelThermalLoadReal

except ImportError:

    class ParallelThermalLoadReal(UnavailableObject):
        pass


try:
    from libaster import ParallelThermalLoadFunction

except ImportError:

    class ParallelThermalLoadFunction(UnavailableObject):
        pass


try:
    from libaster import ParallelDOFNumbering
except ImportError:

    class ParallelDOFNumbering(UnavailableObject):
        pass


try:
    from libaster import ParallelEquationNumbering
except ImportError:

    class ParallelEquationNumbering(UnavailableObject):
        pass


try:
    from libaster import MGISBehaviour
except ImportError:

    class MGISBehaviour(UnavailableObject):
        pass
