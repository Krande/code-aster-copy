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
:py:class:`ParallelMesh` --- Assignment of parallel mesh
************************************************************************
"""
from collections import Counter

from ..Utilities import ExecutionParameter, Options, injector, logger, MPI
from ..Messages import UTMESS
from .datastructure_ext import OnlyParallelObject

try:
    from libaster import ParallelMesh, Mesh

    from ..Utilities.MedUtils.MEDPartitioner import MEDPartitioner

    HAS_MEDCOUPLING = True
except ImportError:
    HAS_MEDCOUPLING = False

    class ParallelMesh(OnlyParallelObject):
        pass


try:
    from libaster import ConnectionMesh
except ImportError:

    class ConnectionMesh(OnlyParallelObject):
        pass


@injector(ParallelMesh)
class ExtendedParallelMesh:
    cata_sdj = "SD.sd_maillage.sd_maillage"
    orig_init = ParallelMesh.__init__

    def __init__(self):
        self.orig_init()
        if not ExecutionParameter().option & Options.HPCMode:
            UTMESS("I", "SUPERVIS_1")
            ExecutionParameter().enable(Options.HPCMode)

    def readMedFile(self, filename, partitioned=False, verbose=0):
        """Read a MED file containing a mesh and eventually partition it

        Arguments:
            filename (string): name of the MED file
            partitioned (bool) : False if the mesh is not yet partitioned and have to
            be partitioned before readinf
            verbose (int) : 0 - warnings
                            1 - informations about main steps
                            2 - informations about all steps

        Returns:
            bool: True if reading and partionning is ok
        """
        if not HAS_MEDCOUPLING:
            logger.info("The MEDCoupling module is required")
            return

        if partitioned:
            filename_partitioned = filename
        else:
            ms = MEDPartitioner(filename)
            ms.partitionMesh(verbose)
            ms.writeMesh()
            filename_partitioned = ms.writedFilename()

        return self._readPartitionedMedFile(filename_partitioned)

    def checkConsistency(self, filename):
        """Check that the partitioned mesh is consistent, i.e. that all nodes,
           cells and groups are present.

           This method is memory consumption since the sequential mesh is loaded

        Arguments:
            filename (string): name of the full MED file

        Returns:
            bool: True if the partitioned mesh is consistent
        """

        # read std mesh
        mesh = Mesh()
        mesh.readMedFile(filename)

        # tests
        group_no_std = mesh.getGroupsOfNodes(local=False)
        group_no_gl  = self.getGroupsOfNodes(local=False)
        if sorted(group_no_std) != sorted(group_no_gl):
            return False

        group_ma_std = mesh.getGroupsOfCells(local=False)
        group_ma_gl  = self.getGroupsOfCells(local=False)
        if sorted(group_ma_std) != sorted(group_ma_gl):
            return False

        nb_nodes_std = mesh.getNumberOfNodes()
        nb_nodes_lc = len(self.getInnerNodes())
        nb_nodes_gl = MPI.COMM_WORLD.allreduce(nb_nodes_lc, MPI.SUM)
        if nb_nodes_std != nb_nodes_gl:
            return False

        rank = MPI.COMM_WORLD.Get_rank()
        nb_cells_std = mesh.getNumberOfCells()
        cells_rank = self.getCellsRank()
        nb_cells_lc = Counter(cells_rank)[rank]
        nb_cells_gl = MPI.COMM_WORLD.allreduce(nb_cells_lc, MPI.SUM)
        if nb_cells_std != nb_cells_gl:
            return False

        return True

@injector(ConnectionMesh)
class ExtendedConnectionMesh:
    cata_sdj = "SD.sd_maillage.sd_connection_mesh"
