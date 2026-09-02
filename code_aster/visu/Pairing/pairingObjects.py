# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
from .meshMatplotlibFigure import meshMatplotlibFigure as MMFig
from abc import ABC, abstractmethod
import itertools

AMBIENT_DIMENSION = [2, 3]


## -----------------------------------------------------------
#   GENERIC CLASS FOR PAIRING OBJECTS
## -----------------------------------------------------------
class PairingObject(ABC):
    def __init__(self, dimension, masterDomain, masterInterface, slaveDomain, slaveInterface):
        """Constructor

        Arguments
        ---------
        dimension : int
            Dimension of the problem (2 or 3).
        masterDomain : str
            Name of the master solid.
        masterInterface : str
            Name of the master interface.
        slaveDomain : str
            Name of the slave solid.
        slaveInterface : str
            Name of the slave interface.
        """
        assert dimension in AMBIENT_DIMENSION
        self._dim = dimension
        # - mesh informations
        self._mastSolid = masterDomain
        self._mastInt = masterInterface
        self._slvSolid = slaveDomain
        self._slvtInt = slaveInterface
        # - List and arrays needed
        self._coords = None
        self._asterConnectivity = None
        self._listPairs = None
        self._listIntersectionPts = None
        self._listQuadraturePts = None
        self._listPairsBasicInfo = None
        self._listPairsDict = None
        # - Flag for functionnalities
        self._flag_MeshInfos = False
        self._flag_PairingInfos = False
        self._flag_CellsInfos = False

    def getNodesCoordsFromCellIndices(self, cellIndices):
        """Provide a list of nodes belonging to a given list of cells

        Arguments
        ---------
        cellIndices : list
            Cells indices
        """
        if self._flag_MeshInfos:
            nodesIndices = [
                self._asterConnectivity[ind] for ind in cellIndices
            ]  # self._asterConnectivity[cellIndices]
            nodesIndices = list(itertools.chain(*nodesIndices))
            nodesIndices = list(set(nodesIndices))
            nodesIndices.sort()

            return self._coords[nodesIndices], nodesIndices
        else:
            raise ValueError("Mesh informations have not been implemented.")

    def checkInfosForPlot(self):
        """Check if all the infos necessary for plots have been provided"""
        boolTest = self._flag_MeshInfos and self._flag_PairingInfos and self._flag_CellsInfos
        if not boolTest:
            raise ValueError(
                f"Pairing Object not properly initialized for plot, Mesh info:{self._flag_MeshInfos}, Pairing Info:{self._flag_PairingInfos}, Cells Info: {self._flag_CellsInfos}"
            )

    def plotInteractive(
        self,
        optionMesh,
        suboptionMesh,
        optionPair,
        addNodeLabel=False,
        addMeshNodes=False,
        addLegend=False,
        index=None,
        s=50,
        indexPlaneProjected=None,
    ):
        """Interative plot

        Arguments
        ---------
        optionMesh : str
            Option for mesh visu (domain or interface).
        suboptionMesh : str
            Target problem (all pairs ? one selected pair ?).
        optionPair : str
            Option for the level of information one seek (mesh ? pairs ? intersection points ? quadrature points ?).
        addNodeLabel : bool
            If True, then add node labels (numbering).
        addMeshNodes : bool
            If True, then add mesh node (bullet for nodes).
        addLegend : bool
            If True, add the legend.
        index : str
            Index of a selected cell.
        s : float
            Opacity parameter.
        indexPlaneProjected : str
            Option for projection when dealing with 3D case.
        """
        self.checkInfosForPlot()
        fig = MMFig(
            self,
            self._dim,
            optionMesh,
            suboptionMesh,
            optionPair,
            addNodeLabel=addNodeLabel,
            addMeshNodes=addMeshNodes,
            addLegend=addLegend,
            index=index,
            indexPlaneProjected=indexPlaneProjected,
        )
        fig.plot(s)

    def computebasicInfosFromPairs(self):
        """Compute information associated to each cell in all the provided pairs"""
        # - Step 1: Compute unique indices for the first column
        unique_indices, counts = np.unique(self._listPairs[:, 0], return_counts=True)
        # Result as list [[index, number of occurrences]]
        result_counts = list(zip(unique_indices, counts))
        # - Step 2
        index_dict = {}

        for i, unique_index in enumerate(unique_indices):
            # Find all pairs where unique_index is present in the first column
            indices_pairs = np.where(self._listPairs[:, 0] == unique_index)[0]
            indices_cell = self._listPairs[indices_pairs, 1].tolist()
            # Fill dictionnary
            index_dict[unique_index] = {
                "indicesCell": indices_cell,
                "indicesPairs": indices_pairs.tolist(),
            }

        # - Save data
        self._listPairsBasicInfo = np.copy(result_counts)
        self._listPairsDict = index_dict

    def getSlaveCellsPaired(self):
        """Return the indices of slave cells that appear in pairs"""
        if self._listPairsBasicInfo is None:
            self.computebasicInfosFromPairs()
        return self._listPairsBasicInfo[:, 0]

    @abstractmethod
    def setMeshInfos(self, *args, **kwargs):
        """Set mesh informations: node coordinates and connectivity"""
        pass

    @abstractmethod
    def setPairingInfos(self, *args, **kwargs):
        """Set pairing informations:
        list of pairs of cells
        list of intersection points
        list of quadrature points"""
        pass

    @abstractmethod
    def setCellInfos(self, *args, **kwargs):
        """Set call informations:
        list of indices for the slave cells (domain and interface)
        list of indices for the slave cells (domain and interface)"""
        pass


## -----------------------------------------------------------
#   DERIVED CLASS DEPENDING OF PREVIOUSLY COMPUTED DATA
## -----------------------------------------------------------
# - Class to use when data dumped from pairing process
class PairingAnalysisAsterFromPkl(PairingObject):
    def __init__(self, dimension, masterDomain, masterInterface, slaveDomain, slaveInterface):
        """Constructor

        Arguments
        ---------
        dimension : int
            Dimension of the problem (2 or 3).
        masterDomain : str
            Name of the master solid.
        masterInterface : str
            Name of the master interface.
        slaveDomain : str
            Name of the slave solid.
        slaveInterface : str
            Name of the slave interface.
        """
        super().__init__(dimension, masterDomain, masterInterface, slaveDomain, slaveInterface)

    def setMeshInfos(self, coords, asterConnectivity):
        """Set mesh informations: node coordinates and connectivity

        Arguments
        ---------
        coords : list
            Node coordinates.
        asterConnectivity : list
            Connectivity array.
        """
        self._coords = np.copy(coords)
        self._asterConnectivity = asterConnectivity
        # - Update flag
        self._flag_MeshInfos = True

    def setPairingInfos(self, listPairs, listIntersectionPts, listQuadraturePts):
        """Set pairing informations:

        Arguments
        ---------
        listPairs : list
            List of pairs of cells.
        listIntersectionPts : list
            List of intersection points.
        listQuadraturePts : list
            List of quadrature points.
        """
        # - Set pairing information
        self._listPairs = np.copy(listPairs)
        self._listIntersectionPts = np.copy(listIntersectionPts)
        self._listQuadraturePts = np.copy(listQuadraturePts)
        # - Update flag
        self._flag_PairingInfos = True

    def setCellInfos(
        self, indices_grma_slv, indices_grma_mas, indices_grma_do_slv, indices_grma_do_mas
    ):
        """Set cell informations:

        Arguments
        ---------
        indices_grma_slv : list
            List of indices for the slave cells (interface).
        indices_grma_mas : list
            List of indices for the master cells (interface).
        indices_grma_do_slv : list
            List of indices for the slave cells (domain).
        indices_grma_do_mas : list
            List of indices for the master cells (domain).
        """
        # - Slave side
        self._indicesSlaveDomain = indices_grma_do_slv
        self._indicesSlaveInterface = indices_grma_slv
        # - Master side
        self._indicesMasterDomain = indices_grma_do_mas
        self._indicesMasterInterface = indices_grma_mas
        # - Update flag
        self._flag_CellsInfos = True


# - Class to use within a monolithic framework
class PairingAnalysisAster(PairingObject):
    def __init__(
        self,
        dimension,
        masterDomain,
        masterInterface,
        slaveDomain,
        slaveInterface,
        asterPairingProcess,
    ):
        """Constructor

        Arguments
        ---------
        dimension : int
            Dimension of the problem (2 or 3).
        masterDomain : str
            Name of the master solid.
        masterInterface : str
            Name of the master interface.
        slaveDomain : str
            Name of the slave solid.
        slaveInterface : str
            Name of the slave interface.
        asterPairingProcess : AsterPairingProcess
            Previously computed AsterPairingProcess object.
        """
        super().__init__(dimension, masterDomain, masterInterface, slaveDomain, slaveInterface)
        self.setMeshInfos(asterPairingProcess)
        self.setPairingInfos(asterPairingProcess)
        self.setCellInfos(asterPairingProcess)

    def setMeshInfos(self, asterPairingProcess):
        """Set mesh informations from an AsterPairingProcess object.
        Mesh information include:
            node coordinates
            connectivity

        Arguments
        ---------
        asterPairingProcess : AsterPairingProcess
            Previously computed AsterPairingProcess object.
        """
        self._coords = np.copy(asterPairingProcess._coords)
        self._asterConnectivity = asterPairingProcess._asterConnectivity
        # - Update flag
        self._flag_MeshInfos = True

    def setPairingInfos(self, asterPairingProcess):
        """Set pairing informations from an AsterPairingProcess object.
        Pairing information include:
            list of pairs of cells
            list of intersection points
            list of quadrature points

        Arguments
        ---------
        asterPairingProcess : AsterPairingProcess
            Previously computed AsterPairingProcess object.
        """
        if asterPairingProcess._hasRun:
            self._listPairs = asterPairingProcess._listPairs
            self._listIntersectionPts = asterPairingProcess._intePointsList
            self._listQuadraturePts = asterPairingProcess._quadPointsList
            # - Update flag
            self._flag_PairingInfos = True
        else:
            raise ValueError("No pairing has been computed before")

    def setCellInfos(self, asterPairingProcess):
        """Set cell informations from an AsterPairingProcess object.
        Cell informations include:
            list of indices for the slave cells (domain and interface)
            list of indices for the master cells (domain and interface)

        Arguments
        ---------
        asterPairingProcess : AsterPairingProcess
            Previously computed AsterPairingProcess object.
        """
        # - Slave side
        self._indicesSlaveDomain = asterPairingProcess._asterMesh.getCells(self._slvSolid)
        self._indicesSlaveInterface = asterPairingProcess._asterMesh.getCells(self._slvtInt)
        # - Master side
        self._indicesMasterDomain = asterPairingProcess._asterMesh.getCells(self._mastSolid)
        self._indicesMasterInterface = asterPairingProcess._asterMesh.getCells(self._mastInt)

        def transformInt(_list):
            return [int(val) for val in _list]

        self._indicesSlaveDomain = transformInt(self._indicesSlaveDomain)
        self._indicesSlaveInterface = transformInt(self._indicesSlaveInterface)
        self._indicesMasterDomain = transformInt(self._indicesMasterDomain)
        self._indicesMasterInterface = transformInt(self._indicesMasterInterface)
        # - Update flag
        self._flag_CellsInfos = True
