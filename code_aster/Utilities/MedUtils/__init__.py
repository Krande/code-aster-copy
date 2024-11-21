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
This module gives common utilities.

No external import of other :py:mod:`code_aster` packages.
"""

from ..base_utils import config

# aslint: disable=C4008
if config.get("ASTER_HAVE_MPI") and config.get("ASTER_HAVE_MED"):
    from ...Objects import MedFileReader, IncompleteMesh, MeshBalancer, MeshConnectionGraph
    from ...Objects import MedFileAccessType
    from ...Objects import PtScotchPartitioner
    from ...Objects import FieldCharacteristics, SimpleFieldOnNodesReal, Result
    from ...Objects import SimpleFieldOnCellsReal
    from ...Objects import ParallelMesh
    from .medtoasterconnectivity import MYMED2ASTER_CONNECT, MED_TYPES, ASTER_TYPES, toAsterGeoType
else:
    MedFileReader = IncompleteMesh = MeshBalancer = MeshConnectionGraph = PtScotchPartitioner = None
    FieldCharacteristics = SimpleFieldOnNodesReal = Result = SimpleFieldOnCellsReal = None
    MYMED2ASTER_CONNECT = MED_TYPES = ASTER_TYPES = MedFileAccessType = toAsterGeoType = None
    ParallelMesh = None
