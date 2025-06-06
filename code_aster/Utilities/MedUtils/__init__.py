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

# person_in_charge: mathieu.courtois@edf.fr

"""
This module gives common utilities.

No external import of other :py:mod:`code_aster` packages.
"""

from ..base_utils import config

# aslint: disable=C4008
MedFileReader = MedFileAccessType = None
MYMED2ASTER_CONNECT = MED_TYPES = ASTER_TYPES = toAsterGeoType = None
PtScotchPartitioner = None
IncompleteMesh = MeshBalancer = MeshConnectionGraph = ParallelMesh = Model = None
FieldCharacteristics = SimpleFieldOnCellsReal = SimpleFieldOnNodesReal = Result = None

if config.get("ASTER_HAVE_MPI"):
    if config.get("ASTER_HAVE_MED"):
        from ...Objects import MedFileAccessType, MedFileReader
        from .medtoasterconnectivity import (
            ASTER_TYPES,
            MED_TYPES,
            MYMED2ASTER_CONNECT,
            toAsterGeoType,
        )
    if config.get("ASTER_HAVE_PTSCOTCH"):
        from ...Objects import PtScotchPartitioner
    from ...Objects import (
        FieldCharacteristics,
        IncompleteMesh,
        MeshBalancer,
        MeshConnectionGraph,
        ParallelMesh,
        Result,
        SimpleFieldOnCellsReal,
        SimpleFieldOnNodesReal,
        Model,
    )
