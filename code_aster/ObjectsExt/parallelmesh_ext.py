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

# person_in_charge: nicolas.pignet@edf.fr
"""
:py:class:`ParallelMesh` --- Assignment of parallel mesh
************************************************************************
"""

import os.path as osp

import numpy as np

from ..Commands import CREA_MAILLAGE
from ..Messages import UTMESS
from ..Objects import (
    CommGraph,
    ConnectionMesh,
    IncompleteMesh,
    Mesh,
    ParallelMesh,
    PythonBool,
    ResultNaming,
)
from ..Objects.Serialization import InternalStateBuilder
from ..Utilities import MPI, ExecutionParameter, Options, force_list, injector, shared_tmpdir
from ..Utilities.MedUtils.MedMeshAndFieldsSplitter import splitMeshAndFieldsFromMedFile
from . import mesh_builder


class ParallelMeshStateBuilder(InternalStateBuilder):
    """Class that returns the internal state of a *ParallelMesh*."""

    def restore(self, mesh):
        """Restore the *DataStructure* content from the previously saved internal
        state.

        Arguments:
            mesh (*DataStructure*): The *DataStructure* object to be restored.
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

    def readMedFile(self, filename, meshname=None, partitioned=False, verbose=0):
        """Read a MED file containing a mesh and eventually partition it.

        Arguments:
            filename (str): Name of the MED file.
            meshname (str): Name of the mesh to be read from file.
            partitioned (bool): False if the mesh is not yet partitioned and have to
                be partitioned before reading.
            verbose (int): Verbosity between 0 (a few details) to 2 (more verbosy).
        """
        if not partitioned:
            self, field = splitMeshAndFieldsFromMedFile(filename, outMesh=self)
            self.show(verbose & 3)
        else:
            mesh_builder.buildFromMedFile(self, filename, meshname, verbose)

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
        rank = MPI.ASTER_COMM_WORLD.Get_rank()

        nb_nodes_lc = len(self.getInnerNodes())
        nb_nodes_gl = MPI.ASTER_COMM_WORLD.allreduce(nb_nodes_lc, MPI.SUM)

        nb_cells_lc = len(self.getInnerCells())
        nb_cells_gl = MPI.ASTER_COMM_WORLD.allreduce(nb_cells_lc, MPI.SUM)

        test = True

        if rank == 0:
            mesh = Mesh()
            mesh.readMedFile(filename)

            # tests
            group_no_std = mesh.getGroupsOfNodes(local=False)
            group_no_gl = self.getGroupsOfNodes(local=False)
            test1 = sorted(group_no_std) == sorted(group_no_gl)
            test = test and test1
            if not test1:
                print(
                    f"FAILED {filename} Groups of nodes differ:",
                    sorted(group_no_std),
                    sorted(group_no_gl),
                )

            group_ma_std = mesh.getGroupsOfCells(local=False)
            group_ma_gl = self.getGroupsOfCells(local=False)
            test1 = sorted(group_ma_std) == sorted(group_ma_gl)
            test = test and test1
            if not test1:
                print(
                    f"FAILED {filename} Groups of cells differ:",
                    sorted(group_ma_std),
                    sorted(group_ma_gl),
                )

            nb_nodes_std = mesh.getNumberOfNodes()
            test1 = nb_nodes_std == nb_nodes_gl
            test = test and test1
            if not test1:
                print(f"FAILED {filename} Number of nodes differs:", nb_nodes_std, nb_nodes_gl)

            nb_cells_std = mesh.getNumberOfCells()
            test1 = nb_cells_std == nb_cells_gl
            test = test and test1
            if not test1:
                print(f"FAILED {filename} Number of cells differs:", nb_cells_std, nb_cells_gl)

        return MPI.ASTER_COMM_WORLD.bcast(test, root=0)

    def checkJoints(self):
        comm = MPI.COMM_WORLD
        rank = MPI.ASTER_COMM_WORLD.Get_rank()
        l2G = self.getLocalToGlobalNodeIds()
        graph = CommGraph()

        dictProc = {}
        cmpt = 0
        for proc in self.getOppositeDomains():
            graph.addCommunication(proc)
            dictProc[proc] = cmpt
            cmpt += 1
        graph.synchronizeOverProcesses()

        outN = self.getOuterNodes()
        dictOutN = {}
        for nId in outN:
            dictOutN[nId + 1] = 1

        coords = self.getCoordinates()

        matchings = graph.getMatchings()
        toReturn = True
        for procT in matchings:
            proc = procT[1]
            tag = procT[0]
            if proc == -1:
                continue
            numJoint = dictProc[proc]
            fJ = self.getSendJoint(numJoint)
            gFJ = []
            curCoordsFJ = []
            for i in range(int(len(fJ) / 2)):
                nId = fJ[2 * i]
                gFJ.append(l2G[nId - 1])
                curCoordsFJ.append(coords[(nId - 1)])

            sJ = self.getReceiveJoint(numJoint)
            gSJ = []
            curCoordsSJ = []
            for i in range(int(len(sJ) / 2)):
                nId = sJ[2 * i]
                dictOutN.pop(nId)
                gSJ.append(l2G[nId - 1])
                curCoordsSJ.append(coords[(nId - 1)])

            if proc < rank:
                comm.send(gFJ, dest=proc, tag=tag)
                data1 = comm.recv(source=proc, tag=tag)
                if data1 != gFJ:
                    print("Rank", rank, "Opposite domain", proc, "Joint number", numJoint, "NOOK")
                    toReturn = False
                comm.send(curCoordsFJ, dest=proc, tag=tag)
                data2 = comm.recv(source=proc, tag=tag)
                tab1 = np.array(curCoordsFJ)
                tab2 = np.array(data2)
                if np.allclose(tab1, tab2, atol=1e-12) is not True:
                    toReturn = False
                    print("Node joint coordinates are not matching")
            else:
                data1 = comm.recv(source=proc, tag=tag)
                comm.send(gSJ, dest=proc, tag=tag)
                if data1 != gSJ:
                    print("Rank", rank, "Opposite domain", proc, "Joint number", numJoint, "NOOK")
                    toReturn = False
                comm.send(curCoordsSJ, dest=proc, tag=tag)
                data2 = comm.recv(source=proc, tag=tag)
                tab1 = np.array(curCoordsSJ)
                tab2 = np.array(data2)
                if np.allclose(tab1, tab2, atol=1e-12) is not True:
                    toReturn = False
                    print("Node joint coordinates are not matching")
        if len(dictOutN) != 0:
            print("Some outer nodes are not in a joint")
            toReturn = False
        return toReturn

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
    def buildSquare(cls, l=1, refine=0, info=1):
        """Build the quadrilateral mesh of a square.

        Arguments:
            l [float] : size of the cube (default 1.).
            refine [int] : number of mesh refinement iterations (default 0).
            info [int] : verbosity mode (1 or 2). (default 1).
        """

        ### Refine some levels on whole mesh, the remaining after partitioning
        min_level = 6
        refine_0 = min(min_level, refine)
        refine_1 = refine - refine_0

        with shared_tmpdir("buildSquare") as tmpdir:
            filename = osp.join(tmpdir, "buildSquare.med")
            if MPI.ASTER_COMM_WORLD.Get_rank() == 0:
                mesh = Mesh.buildSquare(l=l, refine=refine_0, info=info)
                mesh.printMedFile(filename)
            ResultNaming.syncCounter()

            # Mesh creation
            mesh_p = cls()
            mesh_p.readMedFile(filename, verbose=info - 1)
            return mesh_p.refine(refine_1, info)

    @classmethod
    def buildCube(cls, l=1, refine=0, info=1):
        """Build the quadrilateral mesh of a cube.

        Arguments:
            l [float] : size of the cube (default 1.).
            refine [int] : number of mesh refinement iterations (default 0).
            info [int] : verbosity mode (1 or 2). (default 1).
        """

        ### Refine some levels on whole mesh, the remaining after partitioning
        min_level = 7 if MPI.ASTER_COMM_WORLD.Get_size() > 512 else 6
        refine_0 = min(min_level, refine)
        refine_1 = refine - refine_0

        with shared_tmpdir("buildCube") as tmpdir:
            filename = osp.join(tmpdir, "buildCube.med")
            if MPI.ASTER_COMM_WORLD.Get_rank() == 0:
                mesh = Mesh.buildCube(l=l, refine=refine_0, info=info)
                mesh.printMedFile(filename)
            ResultNaming.syncCounter()

            # Mesh creation
            mesh_p = cls()
            mesh_p.readMedFile(filename, verbose=info - 1)
            return mesh_p.refine(refine_1, info)

    @classmethod
    def buildDisk(cls, radius=1, refine=0, info=1):
        """Build the quadrilateral mesh of a disk.

        Arguments:
            radius [float] : radius of the disk (default 1).
            refine [int] : number of mesh refinement iterations (default 0).
            info [int] : verbosity mode (1 or 2). (default 1).
        """

        ### Refine some levels on whole mesh, the remaining after partitioning
        min_level = 6
        refine_0 = min(min_level, refine)
        refine_1 = refine - refine_0

        with shared_tmpdir("buildDisk") as tmpdir:
            filename = osp.join(tmpdir, "buildDisk.med")
            if MPI.ASTER_COMM_WORLD.Get_rank() == 0:
                mesh = Mesh.buildDisk(radius=radius, refine=refine_0, info=info)
                mesh.printMedFile(filename)
            ResultNaming.syncCounter()

            # Mesh creation
            mesh_p = cls()
            mesh_p.readMedFile(filename, verbose=info - 1)
            return mesh_p.refine(refine_1, info)

    @classmethod
    def buildCylinder(cls, height=3, radius=1, refine=0, info=1):
        """Build the hexaedral mesh of a cylinder.

        Arguments:
            height [float] : height of the cylinder along the z axis (default 0).
            radius [float] : radius of the cylinder (default 1).
            refine [int] : number of mesh refinement iterations (default 0).
            info [int] : verbosity mode (1 or 2). (default 1).
        """

        ### Refine some levels on whole mesh, the remaining after partitioning
        min_level = 6
        refine_0 = min(min_level, refine)
        refine_1 = refine - refine_0

        with shared_tmpdir("buildCylinder") as tmpdir:
            filename = osp.join(tmpdir, "buildCylinder.med")
            if MPI.ASTER_COMM_WORLD.Get_rank() == 0:
                mesh = Mesh.buildCylinder(height=height, radius=radius, refine=refine_0, info=info)
                mesh.printMedFile(filename)
            ResultNaming.syncCounter()

            # Mesh creation
            mesh_p = cls()
            mesh_p.readMedFile(filename, verbose=info - 1)
            return mesh_p.refine(refine_1, info)

    def getNodes(self, group_name=[], localNumbering=True, same_rank=None):
        """Return the list of the indexes of the nodes that belong to a group of nodes.

        Arguments:
            group_name (str): Name of groups (default: [] = all nodes).
            localNumbering (bool) : use local or global numbering (default: True)
            same_rank : - None: keep all nodes (default: None)
                        - True: keep the nodes which are owned by the current MPI-rank
                        - False: keep the nodes which are not owned by the current MPI-rank

        Returns:
            list[int]: Indexes of the nodes of groups.
        """

        val = {None: PythonBool.NONE, True: PythonBool.TRUE, False: PythonBool.FALSE}

        return self._getNodes(force_list(group_name), localNumbering, val[same_rank])

    def getNodesFromCells(self, group_name, localNumbering=True, same_rank=None):
        """Returns the nodes indexes of a group of cells.

        Arguments:
            group_name (str/list[str]): Name of the group.
            localNumbering (bool) : use local or global numbering (default: True)
            same_rank : - None: keep all nodes (default: None)
                        - True keep the nodes which are owned by the current MPI-rank
                        - False: keep the nodes which are not owned by the current MPI-rank

        Returns:
            list[int]: Indexes of the nodes of the group.
        """

        val = {None: PythonBool.NONE, True: PythonBool.TRUE, False: PythonBool.FALSE}

        return self._getNodesFromCells(force_list(group_name), localNumbering, val[same_rank])


@injector(ConnectionMesh)
class ExtendedConnectionMesh:
    cata_sdj = "SD.sd_maillage.sd_connection_mesh"


@injector(IncompleteMesh)
class ExtendedIncompleteMesh:
    cata_sdj = "SD.sd_maillage.sd_maillage"

    def readMedFile(self, filename, meshname=None, verbose=1):
        """Read a MED file containing a mesh.

        Arguments:
            filename (string): name of the MED file
            meshname (str): Name of the mesh to be read from file.
            verbose (int) : 0 - warnings
                            1 - informations about main steps
                            2 - informations about all steps
        """
        mesh_builder.buildFromMedFile(self, filename, meshname, verbose)
