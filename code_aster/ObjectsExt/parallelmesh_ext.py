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

# person_in_charge: nicolas.pignet@edf.fr
"""
:py:class:`ParallelMesh` --- Assignment of parallel mesh
************************************************************************
"""
import os.path as osp

from ..Commands import CREA_MAILLAGE
from ..Messages import UTMESS
from ..Objects import ConnectionMesh, Mesh, ParallelMesh, ResultNaming, PythonBool
from ..Objects.Serialization import InternalStateBuilder
from ..Utilities import MPI, ExecutionParameter, Options, injector, logger, shared_tmpdir

try:
    from ..Utilities.MedUtils.MEDPartitioner import MEDPartitioner

    HAS_MEDCOUPLING = True
except ImportError:
    HAS_MEDCOUPLING = False


class ParallelMeshStateBuilder(InternalStateBuilder):
    """Class that returns the internal state of a *ParallelMesh*."""

    def restore(self, mesh):
        """Restore the *DataStructure* content from the previously saved internal
        state.

        Arguments:
            mesh (*DataStructure*): The *DataStructure* object to be pickled.
        """
        super().restore(mesh)
        mesh.build()
        mesh._updateGlobalGroupOfNodes()
        mesh._updateGlobalGroupOfCells()


@injector(ParallelMesh)
class ExtendedParallelMesh:
    cata_sdj = "SD.sd_maillage.sd_maillage"
    internalStateBuilder = ParallelMeshStateBuilder
    orig_init = ParallelMesh.__init__

    def __init__(self, *args, **kwargs):
        self.orig_init(*args, **kwargs)
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
        rank = MPI.COMM_WORLD.Get_rank()

        nb_nodes_lc = len(self.getInnerNodes())
        nb_nodes_gl = MPI.COMM_WORLD.allreduce(nb_nodes_lc, MPI.SUM)

        nb_cells_lc = len(self.getInnerCells())
        nb_cells_gl = MPI.COMM_WORLD.allreduce(nb_cells_lc, MPI.SUM)

        test = True

        if rank == 0:
            mesh = Mesh()
            mesh.readMedFile(filename)

            # tests
            group_no_std = mesh.getGroupsOfNodes(local=False)
            group_no_gl = self.getGroupsOfNodes(local=False)
            test = sorted(group_no_std) == sorted(group_no_gl)

            group_ma_std = mesh.getGroupsOfCells(local=False)
            group_ma_gl = self.getGroupsOfCells(local=False)
            test = test and sorted(group_ma_std) == sorted(group_ma_gl)

            nb_nodes_std = mesh.getNumberOfNodes()
            test = test and nb_nodes_std == nb_nodes_gl

            nb_cells_std = mesh.getNumberOfCells()
            test = test and nb_cells_std == nb_cells_gl

        return MPI.COMM_WORLD.bcast(test, root=0)

    def refine(self, ntimes=1, info=1):
        """Refine the mesh uniformly. Each edge is split in two.

        Arguments:
            ntimes [int] : the number of times the mesh is to be refined.
            info [int] : verbosity mode (1 or 2). Default 1.

        Returns:
            ParallelMesh: the refined mesh.
        """

        return CREA_MAILLAGE(MAILLAGE=self, RAFFINEMENT=_F(TOUT="OUI", NIVEAU=ntimes), INFO=info)

    @classmethod
    def buildSquare(self, lx=1.0, ly=1.0, refine=0, info=1):
        """Build the quadrilateral mesh of a square.

        Arguments:
            lx [float] : length of the square along the x axis (default 1.).
            ly [float] : length of the square along the y axis (default 1.).
            refine [int] : number of mesh refinement iterations (default 0).
            info [int] : verbosity mode (1 or 2). (default 1).
        """

        with shared_tmpdir("buildSquare") as tmpdir:
            filename = osp.join(tmpdir, "buildSquare.med")
            # nrefine is added to create a mesh with enougth cells
            # to be partitioned equally and not generated a too big file
            nrefine = min(9, refine)
            if MPI.COMM_WORLD.Get_rank() == 0:
                mesh = Mesh.buildSquare(
                    lx=lx, ly=ly, refine=nrefine, info=info)
                mesh.printMedFile(filename)
            ResultNaming.syncCounter()

            # Mesh creation
            mesh_p = ParallelMesh()
            mesh_p.readMedFile(filename, verbose=info - 1)

            # Mesh refinement
            nrefinep = refine - nrefine
            return CREA_MAILLAGE(
                MAILLAGE=mesh_p, RAFFINEMENT=_F(TOUT="OUI", NIVEAU=nrefinep), INFO=info
            )

    @classmethod
    def buildCube(self, lx=1.0, ly=1.0, lz=1.0, refine=0, info=1):
        """Build the hexaedral mesh of a cube.

        Arguments:
            lx [float] : length of the cube along the x axis (default 1.).
            ly [float] : length of the cube along the y axis (default 1.).
            lz [float] : length of the cube along the z axis (default 1.).
            refine [int] : number of mesh refinement iterations (default 0).
            info [int] : verbosity mode (1 or 2). (default 1).
        """

        with shared_tmpdir("buildCube") as tmpdir:
            filename = osp.join(tmpdir, "buildCube.med")
            # nrefine is added to create a mesh with enougth cells
            # to be partitioned equally and not generated a too big file
            min_level = 6
            if MPI.COMM_WORLD.Get_size() > 512:
                min_level = 7
            nrefine = min(min_level, refine)
            if MPI.COMM_WORLD.Get_rank() == 0:
                mesh = Mesh.buildCube(
                    lx=lx, ly=ly, lz=lz, refine=nrefine, info=info)
                mesh.printMedFile(filename)
            ResultNaming.syncCounter()

            # Mesh creation
            mesh_p = ParallelMesh()
            mesh_p.readMedFile(filename, verbose=info - 1)

            # Mesh refinement
            nrefinep = refine - nrefine
            return CREA_MAILLAGE(
                MAILLAGE=mesh_p, RAFFINEMENT=_F(TOUT="OUI", NIVEAU=nrefinep), INFO=info
            )

    @classmethod
    def buildDisk(self, radius=1, refine=0, info=1):
        """Build the quadrilateral mesh of a disk.

        Arguments:
            radius [float] : radius of the disk (default 1).
            refine [int] : number of mesh refinement iterations (default 0).
            info [int] : verbosity mode (1 or 2). (default 1).
        """

        with shared_tmpdir("buildDisk") as tmpdir:
            filename = osp.join(tmpdir, "buildDisk.med")
            if MPI.COMM_WORLD.Get_rank() == 0:
                mesh = Mesh.buildDisk(radius, refine, info)
                mesh.printMedFile(filename)
            ResultNaming.syncCounter()

            # Mesh creation
            mesh_p = ParallelMesh()
            mesh_p.readMedFile(filename, verbose=info - 1)
            return mesh_p

    @classmethod
    def buildCylinder(self, height=3, radius=1, refine=0, info=1):
        """Build the hexaedral mesh of a cylinder.

        Arguments:
            height [float] : height of the cylinder along the z axis (default 0).
            radius [float] : radius of the cylinder (default 1).
            refine [int] : number of mesh refinement iterations (default 0).
            info [int] : verbosity mode (1 or 2). (default 1).
        """

        with shared_tmpdir("buildCylinder") as tmpdir:
            filename = osp.join(tmpdir, "buildCylinder.med")
            if MPI.COMM_WORLD.Get_rank() == 0:
                mesh = Mesh.buildCylinder(height, radius, refine, info)
                mesh.printMedFile(filename)
            ResultNaming.syncCounter()

            # Mesh creation
            mesh_p = ParallelMesh()
            mesh_p.readMedFile(filename, verbose=info - 1)
            return mesh_p

    def getNodes(self, group_name="", localNumbering=True, same_rank=None):
        """ Return the list of the indexes of the nodes that belong to a group of nodes.

            Arguments:
                group_name (str): Name of the group (default: "" = all nodes).
                localNumbering (bool) : use local or global numbering (default: True)
                same_rank : - None: keep all nodes (default: None)
                            - True: keep the nodes which are owned by the current MPI-rank
                            - False: keep the nodes which are not owned by the current MPI-rank

            Returns:
                list[int]: Indexes of the nodes of the group.
        """

        val = {None: PythonBool.NONE,
               True: PythonBool.TRUE,
               False: PythonBool.FALSE}

        return self._getNodes(group_name, localNumbering, val[same_rank])

    def getNodesFromCells(self, group_name, localNumbering=True, same_rank=None):
        """ Returns the nodes indexes of a group of cells.

            Arguments:
                group_name (str): Name of the group.
                localNumbering (bool) : use local or global numbering (default: True)
                same_rank : - None: keep all nodes (default: None)
                            - True keep the nodes which are owned by the current MPI-rank
                            - False: keep the nodes which are not owned by the current MPI-rank

            Returns:
                list[int]: Indexes of the nodes of the group.
        """

        val = {None: PythonBool.NONE,
               True: PythonBool.TRUE,
               False: PythonBool.FALSE}

        return self._getNodesFromCells(group_name, localNumbering, val[same_rank])


@injector(ConnectionMesh)
class ExtendedConnectionMesh:
    cata_sdj = "SD.sd_maillage.sd_connection_mesh"
