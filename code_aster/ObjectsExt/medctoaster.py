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

import numpy as np

from ..Messages import UTMESS
from ..Utilities import MPI, Timer, medcoupling as medc

ASTER_TYPES = [1, 2, 4, 6, 7, 9, 11, 12, 14, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]

MED_TYPES = [
    1,
    102,
    103,
    104,
    203,
    206,
    207,
    204,
    208,
    209,
    304,
    310,
    306,
    315,
    318,
    305,
    313,
    308,
    320,
    327,
]

MED2ASTER_CONNECT = {
    "POINT1": [0],
    "SEG2": range(2),
    "TRI3": range(3),
    "QUAD4": range(4),
    "TETRA4": [0, 2, 1, 3],
    "HEXA8": [0, 3, 2, 1, 4, 7, 6, 5],
    "PYRA5": [0, 3, 2, 1, 4],
    "PENTA6": [0, 2, 1, 3, 5, 4],
    "SEG3": range(3),
    "TRI6": range(6),
    "QUAD8": range(8),
    "TETRA10": [0, 2, 1, 3, 6, 5, 4, 7, 9, 8],
    "HEXA20": [0, 3, 2, 1, 4, 7, 6, 5, 11, 10, 9, 8, 16, 19, 18, 17, 15, 14, 13, 12],
    "PYRA13": [0, 3, 2, 1, 4, 8, 7, 6, 5, 9, 12, 11, 10],
    "PENTA15": [0, 2, 1, 3, 5, 4, 8, 7, 6, 12, 14, 13, 11, 10, 9],
    "SEG4": range(4),
    "TRI7": range(7),
    "QUAD9": range(9),
    "PENTA18": [0, 2, 1, 3, 5, 4, 8, 7, 6, 12, 14, 13, 11, 10, 9, 17, 16, 15],
    "HEXA27": [
        0,
        3,
        2,
        1,
        4,
        7,
        6,
        5,
        11,
        10,
        9,
        8,
        16,
        19,
        18,
        17,
        15,
        14,
        13,
        12,
        20,
        24,
        23,
        22,
        21,
        25,
        26,
    ],
}


# Fonction manquante dans medcoupling
def ConvertToMEDFileGeoType(medcoupling_type):
    """Convert a medcoupling mesh type to a medfile type.

    Arguments:
        medcoupling_type (int): medcoupling type.

    Returns:
        int: medfile type.
    """
    if not hasattr(ConvertToMEDFileGeoType, "cache_dict"):
        ConvertToMEDFileGeoType.cache_dict = {
            medc.MEDFileUMesh.ConvertFromMEDFileGeoType(i): i for i in MED_TYPES
        }
    return ConvertToMEDFileGeoType.cache_dict[medcoupling_type]


def toAsterGeoType(medfile_type):
    """Convert a med mesh type to a code_aster type.

    Arguments:
        medfile_type (int): med type.

    Returns:
        int: code_aster type (index in '&CATA.TM.NOMTM').
    """
    if not hasattr(toAsterGeoType, "cache_dict"):
        toAsterGeoType.cache_dict = dict(zip(MED_TYPES, ASTER_TYPES))
    return toAsterGeoType.cache_dict[medfile_type]


class MEDCouplingMeshHelper:
    """This object extracts properties from a medcoupling mesh as expected
    to create or manipulate a code_aster mesh.
    """

    _medc2aster_connect = {}

    def __init__(self):
        self._mesh = None
        self._mesh_name = ""
        self._file_name = ""
        self._dim = 0
        self._coords = None
        self._groups_of_cells = {}
        self._groups_of_nodes = {}
        self._connectivity_aster = []
        self._types = []
        # for the joints
        self._domains = []
        self._djoints = {}
        self._nodesOwner = []
        self._gNumbering = []
        self._timer = Timer(title="medcoupling conversion")

    def setMedCouplingMesh(self, medmesh):
        """Define the medcoupling mesh to be converted.

        Arguments:
            medmesh (medc.MEDFileUMesh): medcoupling mesh object.
        """
        self._mesh = medmesh
        self._mesh_name = medmesh.getName()

    def readMedFile(self, path: str, meshname: str = None):
        """Define the medcoupling mesh by reading a MED file.

        Arguments:
            path (str): Path to the med file.
            meshname (str): Name of the mesh to be read from file.
        """
        self.setMedCouplingMesh(medc.MEDFileUMesh(path, meshname))
        self._file_name = path

    def buildMesh(self, mesh, verbose=0):
        """Convert the MEDCoupling mesh.

        Arguments:
            mesh (Mesh|ParallelMesh): The Mesh object to be filled.
            verbose (int): Verbosity between 0 (a few details) to 2 (more verbosy).

        Returns:
            Mesh: The same Mesh/ParallelMesh object completed.
        """
        assert mesh.getDimension() == 0, "the mesh object is not empty"
        UTMESS("I", "MED_10", valk=self.name)
        timer = self._timer
        with timer("medcoupling parsing"):
            if not mesh.isIncomplete():
                self.parse()
                with timer("extracting joints"):
                    self.extract_joints()
            else:
                self.parseIncomplete()
        if verbose & 4:
            self.debugPrint()
        with timer("nodes, connectivity"):
            mesh._initDefinition(
                max(self._dim, 2),
                self.coordinates,
                self.connectivity,
                self.types,
                len(self._groups_of_cells),
                len(self._groups_of_nodes),
            )
        groups = self.groupsOfCells
        if groups:
            with timer("groups of cells"):
                mesh._addGroupsOfCells(*groups)
        groups = self.groupsOfNodes
        if groups:
            with timer("groups of nodes"):
                mesh._addGroupsOfNodes(*groups)
        if mesh.isParallel():
            with timer("creating joints"):
                mesh._create_joints(
                    self.domains, self.globalNumbering, self.nodesOwner, self.getAllJoints()
                )
        with timer("completion"):
            if not mesh.isIncomplete():
                mesh._endDefinition()
        # cheat code for debugging and detailed time informations: 1989
        with timer("show details"):
            mesh.show(verbose & 3)
        if verbose & 4:
            print(repr(timer), flush=True)
        return mesh

    def debugPrint(self):
        """Print internal content for debugging."""

        def _info_grp(groups, label):
            if not groups:
                return
            names, groups = groups
            cumsize = sum(map(len, groups))
            minsize = min(map(len, groups))
            meansize = cumsize / len(groups)
            maxsize = max(map(len, groups))
            etc = "..." if len(groups) > 10 else ""
            print(f"+ number of groups of {label}: {len(groups)}, first 10: {names[:10]}{etc}")
            print(f"+ cumsize of groups of {label}: {cumsize}")
            print(f"+ size of groups of {label} (min/mean/max): {minsize}/{meansize:.1f}/{maxsize}")

        print(f"+ mesh name: {self.name}")
        print(f"+ number of cells: {self.numberOfCells}")
        print(f"+ number of nodes: {self.numberOfNodes}")
        _info_grp(self.groupsOfCells, "cells")
        _info_grp(self.groupsOfNodes, "nodes")
        print(f"+ domains: {self.domains}")
        for key, jts in self._djoints.items():
            print(f"+ joints: {key}, size {len(jts)}: {jts[:10]}...")
        print(f"+ number of joints: {len(self._djoints) // 2}", flush=True)

    def parse(self):
        """Walk the medcoupling mesh to extract informations."""
        mesh = self._mesh
        self._dim = mesh.getSpaceDimension()

        nbNodes = mesh.getNumberOfNodes()
        nodes_mc = mesh.getCoords().toNumPyArray()
        nodes_mc.shape = (nbNodes, self._dim)
        self._coords = np.zeros(shape=(nbNodes, 3))
        self._coords[:, : self._dim] = nodes_mc[:, : self._dim]

        non_empty_levels = mesh.getNonEmptyLevels()
        level_shift = 0  # Variable pour la creation d'une numérotation globale

        # Sortie groupes
        # Groupes de cells
        self._groups_of_cells = {}
        self._groups_of_nodes = {}

        self._connectivity_aster = []
        self._types = []
        # groups = mesh.getGroupsNames()
        # NB: `mesh.getGroupsOnSpecifiedLev()` is necessary to avoid segfault
        #     without group of nodes/cells

        # Loop sur les niveaux
        for level in sorted(non_empty_levels):
            mesh_level = mesh[level]

            # Loop sur toutes les types du niveau
            types_at_level = mesh_level.getAllGeoTypesSorted()
            for medcoupling_cell_type in types_at_level:

                # Cells du meme type
                cells_current_type = mesh_level.giveCellsWithType(medcoupling_cell_type)

                # Maillage 1 Single Geo Type (connectivité simple)
                mesh_current_type = medc.MEDCoupling1SGTUMesh(mesh_level[cells_current_type])

                # Type MED
                med_current_type = ConvertToMEDFileGeoType(medcoupling_cell_type)

                # Nombre de noeuds du type
                number_of_nodes_current_type = (
                    medc.MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(medcoupling_cell_type)
                )

                # Remise en forme en tableau
                connectivity_current_type = (
                    mesh_current_type.getNodalConnectivity()
                    .toNumPyArray()
                    .reshape(mesh_current_type.getNumberOfCells(), number_of_nodes_current_type)
                )

                if medcoupling_cell_type != medc.NORM_POINT1:
                    connectivity_current_type = connectivity_current_type[
                        :, self.getConnectivityMedToAster(medcoupling_cell_type)
                    ]

                # Shift de 1
                connectivity_current_type += 1

                # Sauvegarde de la connectivité aster au km
                self._connectivity_aster.extend(connectivity_current_type)
                # Sauvegarde du type med
                self._types.extend([med_current_type] * mesh_current_type.getNumberOfCells())

            # Groupes de cells au niveau
            with self._timer(". store groups of cells"):
                groups = mesh.getGroupsOnSpecifiedLev(level)
                for group in groups:
                    if len(group) > 24:
                        UTMESS("A", "MED_7")
                        continue
                    group_cells = mesh.getGroupArr(level, group).deepCopy()
                    if not group_cells:
                        continue
                    group_cells += 1
                    group_cells += level_shift
                    self._groups_of_cells.setdefault(group, [])
                    self._groups_of_cells[group].extend(group_cells.toNumPyArray())

            # Shift pour numeroration globale
            level_shift += mesh_level.getNumberOfCells()

        # Groupes de noeuds
        with self._timer(". store groups of nodes"):
            groups = mesh.getGroupsOnSpecifiedLev(1)
            for group in groups:
                if len(group) > 24:
                    UTMESS("A", "MED_7")
                    continue
                group_nodes = mesh.getGroupArr(1, group).deepCopy()
                if not group_nodes:
                    continue
                group_nodes += 1
                self._groups_of_nodes.setdefault(group, [])
                self._groups_of_nodes[group].extend(group_nodes.toNumPyArray())

    def parseIncomplete(self):
        """Walk the medcoupling mesh to extract informations."""

        filename = self._file_name
        meshName = self._mesh_name

        rank = MPI.ASTER_COMM_WORLD.Get_rank()
        size = MPI.ASTER_COMM_WORLD.Get_size()
        ttps = medc.GetUMeshGlobalInfo(filename, meshName)[0]
        nodeNb = medc.GetUMeshGlobalInfo(filename, meshName)[3]
        nbNodesProc = int(nodeNb / size)
        startNode = nbNodesProc * rank
        endNode = nbNodesProc * (rank + 1)
        if rank == size - 1:
            endNode = nodeNb

        coords = medc.MEDFileUMesh.LoadPartCoords(
            filename, meshName, -1, -1, ["X", "Y", "Z"], startNode, endNode
        )

        nbNodes = endNode - startNode
        nodes_mc = coords[0].toNumPyArray()
        nodes_mc.shape = (nbNodes, 3)
        self._coords = np.zeros(shape=(nbNodes, 3))
        self._coords[:, :3] = nodes_mc[:, :3]

        params = []
        cts = []
        for tps in ttps:
            for cellType, nbCellsType in tps:
                slc = medc.DataArray.GetSlice(slice(0, nbCellsType, 1), rank, size)
                params += [slc.start, slc.stop, slc.step]
                cts.append(cellType)
                pass
        mrs = medc.MEDFileMeshReadSelector()
        mrs.setNumberOfCoordsLoadSessions(10)
        medFileUMesh = medc.MEDFileUMesh.LoadPartOf(filename, meshName, cts, params, -1, -1, mrs)

        non_empty_levels = medFileUMesh.getNonEmptyLevels()
        level_shift = 0  # Variable pour la creation d'une numérotation globale

        # Sortie groupes
        # Groupes de cells
        self._groups_of_cells = {}
        self._groups_of_nodes = {}

        self._connectivity_aster = []
        self._types = []

        # Loop sur les niveaux
        for level in sorted(non_empty_levels):
            mesh_level = medFileUMesh[level]
            groups = medFileUMesh.getGroupsOnSpecifiedLev(level)
            for group in groups:
                if len(group) > 24:
                    UTMESS("A", "MED_7")
                    continue
                group_cells = medFileUMesh.getGroupArr(level, group)
                group_cells += 1
                group_cells += level_shift
                self._groups_of_cells.setdefault(group, [])
                self._groups_of_cells[group].extend(group_cells.toNumPyArray())

            # Loop sur toutes les types du niveau
            types_at_level = mesh_level.getAllGeoTypesSorted()
            for medcoupling_cell_type in types_at_level:

                # Cells du meme type
                cells_current_type = mesh_level.giveCellsWithType(medcoupling_cell_type)

                # Maillage 1 Single Geo Type (connectivité simple)
                mesh_current_type = medc.MEDCoupling1SGTUMesh(mesh_level[cells_current_type])

                # Type MED
                med_current_type = ConvertToMEDFileGeoType(medcoupling_cell_type)

                # Nombre de noeuds du type
                number_of_nodes_current_type = (
                    medc.MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(medcoupling_cell_type)
                )

                # Remise en forme en tableau
                connectivity_current_type = (
                    mesh_current_type.getNodalConnectivity()
                    .toNumPyArray()
                    .reshape(mesh_current_type.getNumberOfCells(), number_of_nodes_current_type)
                )

                # if medcoupling_cell_type != medc.NORM_POINT1:
                # connectivity_current_type = connectivity_current_type[
                #:, MEDCouplingMeshHelper.getConnectivityMedToAster(medcoupling_cell_type)
                # ]

                # Shift de 1
                connectivity_current_type += 1

                # Sauvegarde de la connectivité aster au km
                self._connectivity_aster.extend(connectivity_current_type)
                self._types.extend([med_current_type] * mesh_current_type.getNumberOfCells())
            level_shift += mesh_level.getNumberOfCells()

        groups = medFileUMesh.getGroupsOnSpecifiedLev(1)
        for group in groups:
            if len(group) > 24:
                UTMESS("A", "MED_7")
                continue
            group_nodes = medFileUMesh.getGroupArr(1, group).deepCopy()
            if not group_nodes:
                continue
            group_nodes += 1
            self._groups_of_nodes.setdefault(group, [])
            self._groups_of_nodes[group].extend(group_nodes.toNumPyArray())

    def extract_joints(self):
        """Extract informations about joints between domains."""
        # _testrank is reserved for testcases, only with one process
        rank = MPI.ASTER_COMM_WORLD.Get_rank()
        self._nodesOwner = np.ones(self._mesh.getNumberOfNodes(), int) * rank
        domains = set()
        self._djoints = {}
        joints = self._mesh.getJoints()
        if joints:
            for name in joints.getJointsNames():
                joint = joints.getJointWithName(name)
                remote = joint.getDomainNumber()
                # assert remote < size and remote != rank # only in parallel!
                domains.add(remote)
                d_in, d_out = [int(i) for i in name.split()]
                assert rank in (
                    d_in,
                    d_out,
                ), f"inconsistent domains: {rank} not in ({d_in}, {d_out})"
                key = ("R" if d_in == rank else "E") + f"{remote:x}".upper()
                assert joint.getNumberOfSteps() == 1, "unexpected value"
                assert joint[0].getNumberOfCorrespondences() == 1, "unexpected value"
                values = joint[0][0].getCorrespondence().toNumPyArray().ravel()
                self._djoints[key] = values.tolist()
                if d_in == rank:
                    values.shape = (len(values) // 2, 2)
                    self._nodesOwner[values[:, 0] - 1] = -1

        self._domains = sorted(domains)
        globNumb = self._mesh.getGlobalNumFieldAtLevel(1)
        if globNumb:
            self._gNumbering = globNumb.toNumPyArray().tolist()

    def getAllJoints(self):
        """Return the definition of all the joints.

        Returns:
            list[list[int]]: Contains the definition of the joints E/R for the
            first domain, E/R for the second, etc. (`2 * number of domains`
            lists of integers).
        """
        if not self._domains:
            return []
        retlist = []
        for dom in self._domains:
            for typ in "ER":
                key = f"{typ}{dom:x}".upper()
                retlist.append(self._djoints[key])
        return retlist

    @property
    def domains(self):
        """list[int]: List of the remote domains."""
        return self._domains

    @property
    def globalNumbering(self):
        """list[int]: List of the global number of each node."""
        return self._gNumbering

    @property
    def nodesOwner(self):
        """list[int]: List of the owner of each node."""
        return self._nodesOwner

    @property
    def mesh(self):
        """medc.MEDFileUMesh: The medcoupling mesh."""
        return self._mesh

    @property
    def name(self):
        """str: The mesh name."""
        return self._mesh_name

    @property
    def numberOfNodes(self):
        """int: Number of nodes."""
        return self._coords.shape[0]

    @property
    def numberOfCells(self):
        """int: Number of cells."""
        return len(self._connectivity_aster)

    @property
    def coordinates(self):
        """list[float]: List of coordinates of each node."""
        return self._coords.flatten().tolist()

    @property
    def connectivity(self):
        """list[list[int]]: Connectivity as a list of lists."""
        return [cell.tolist() for cell in self._connectivity_aster]

    @property
    def medtypes(self):
        """list[int]: List of med type for each cell."""
        return self._types

    @property
    def types(self):
        """list[int]: List of code_aster type for each cell."""
        return [toAsterGeoType(i) for i in self._types]

    @property
    def groupsOfCells(self):
        """dict: Dict of cells contained in each group of cells."""
        return list(zip(*self._groups_of_cells.items()))

    @property
    def groupsOfNodes(self):
        """dict: Dict of nodes contained in each group of nodes."""
        return list(zip(*self._groups_of_nodes.items()))

    @classmethod
    def getConnectivityMedToAster(cls, medcoupling_type):
        """Return the code_aster connectivity of a medcoupling mesh type.

        Arguments:
            medcoupling_type (int): medcoupling type.

        Returns:
            list[int]: List of nodes.
        """
        if not cls._medc2aster_connect:
            cls._medc2aster_connect = {
                getattr(medc, "NORM_%s" % elem): nodes for elem, nodes in MED2ASTER_CONNECT.items()
            }
        return cls._medc2aster_connect[medcoupling_type]
