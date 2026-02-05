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
from ...Objects import MeshPairing
import os
import numpy as np
import pickle
from libaster import (
    PairingMethod,
)  # ContactPairing, ContactComputation, PairingMethod, CoordinatesSpace

## -----------------------------------------------------------
#   AVAILABLE METHODS FOR PAIRING AND MORTAR COMPUTATIONS
## -----------------------------------------------------------
availablePairingMethod = ["BrutForce", "Fast", "Legacy"]


class AsterPairingProcess:

    def __init__(self, groupMaSlv, groupMaMas, asterMesh):
        r"""Constructor

        Args:
            groupMaSlv (:class:`str`): Name of the group of the contact slave interface.
            groupMaMas (:class:`str`): Name of the group of the contact master interface.
            asterMesh (:class:`libaster.Mesh`): aster Mesh considered for pairing

        """
        self._groupMaSlv = groupMaSlv
        self._groupMaMas = groupMaMas
        self._asterMesh = asterMesh
        self._method = None
        # - datastructures initialize to None
        self._coords = None
        self._asterConnectivity = None
        self._listPairs = None
        self._intePointsList = None
        self._quadPointsList = None
        self._hasRun = False

    def setMethod(self, method="BrutForce"):
        r"""Choose the pairing method. Mandatory to run the computation.

        Args:
            method (:class:`str`): Name of the pairing method used in `code_aster`.
            Should be chosen between BrutForce, Legacy, and Fast.

        """
        if method in availablePairingMethod:
            self._method = method
            self._hasRun = False  # the method has changed, then no pairing has been computed
        else:
            raise NameError("method not available: choose between Fast-Legacy-BrutForce")

    def computePairing(self):
        r"""Compute pairing procedure"""

        if self._method is None:
            raise NameError("No pairing method is set: use setMethod before")
        else:

            # - Set MeshPairing and compute pairing
            meshPair = MeshPairing()
            meshPair.setMesh(self._asterMesh)
            meshPair.setPair(self._groupMaSlv, self._groupMaMas)
            if self._method == "BrutForce":
                meshPair.setMethod(PairingMethod.BrutForce)
            elif self._method == "Fast":
                meshPair.setMethod(PairingMethod.Fast)
            elif self._method == "Legacy":
                meshPair.setMethod(PairingMethod.Legacy)

            meshPair.compute()
            # - Save MeshPairing
            self._meshPair = meshPair

    def extractData(self):
        r"""Extract the list of pairs, the list of intersection points and the computed quadrature points coordinates.
        Store them within the datastructure"""
        # - Get the list of pairs
        listPairs = self._meshPair.getListOfPairs()
        # - Get the intersection points
        nbPairs = self._meshPair.getNumberOfPairs()
        intePointsList = []
        for iPair in range(nbPairs):
            IntePoints = self._meshPair.getIntersectionPoints(iPair)
            intePts_current = [tuple(intePt) for intePt in IntePoints]
            intePointsList.append(intePts_current)
        # - Get the quadrature points
        quadPointsList = []
        for iPair in range(nbPairs):
            quadPoints = self._meshPair.getQuadraturePoints(iPair)
            quadPts_current = [tuple(quaPt) for quaPt in quadPoints]
            quadPointsList.append(quadPts_current)
        # Save all data
        self._listPairs = listPairs
        self._intePointsList = intePointsList
        self._quadPointsList = quadPointsList

    def run(self):
        r"""Main method: compute the pairing and save the needed data"""
        # - Compute Pairing
        self.computePairing()
        # - Extract available data
        self.extractData()
        # - Flag to show that one pairing process has been computed
        self._hasRun = True

    def extractMeshInfosFromAsterMesh(self):
        r"""Extract the connecitivity and the node coordinates from the mesh"""
        # - get aster connectivity
        asterConnectivity = self._asterMesh.getConnectivity()
        nbNodes = self._asterMesh.getNumberOfNodes()
        # - get node coordinates for each cell
        coordsFull = self._asterMesh.getCoordinates().getValues()
        coords = np.reshape(coordsFull, (nbNodes, 3))
        # - Save data
        self._coords = coords
        self._asterConnectivity = asterConnectivity

    def computeIndicesMeshFromGroup(self, groupName):
        r"""Returns the indices of the cells belonging to a given group
        Args:
            groupName (:class:`str`): name of the group
        Returns:
            indices (:class:`numpy.ndarray`): list of the cells indices in the group
        """
        return self._asterMesh.getCells(groupName)

    def dumpDataGroupCells(self, groupName, repoSave, nameFile):
        r"""Dump the cell indices of a group in a given directory
        Args:
            groupName (:class:`str`): name of the group
            repoSave (:class:`str`): name of the directory in which we want to save the file
            nameFile (:class:`str`): name of the saved file
        """
        indices = self.computeIndicesMeshFromGroup(groupName)
        if not os.path.exists(repoSave):
            os.makedirs(repoSave)
        with open(os.path.join(repoSave, nameFile + ".pkl"), "wb") as f:
            pickle.dump(indices, f)

    def dumpData(self, repoSave, pairData=True, meshData=True):
        r"""Dump the pairing data in an given directory
        Args:
            groupName (:class:`str`): name of the directory in which we want to save the file
            groupName (:class:`bool`): if True, save pairing Data
            groupName (:class:`str`): if True, save mesh Data (connectivity and node coordinates)
        """
        if not os.path.exists(repoSave):
            os.makedirs(repoSave)
        if pairData:
            if self._coords is not None:
                with open(os.path.join(repoSave, "coords.pkl"), "wb") as f:
                    pickle.dump(self._coords, f)
            if self._asterConnectivity is not None:
                with open(os.path.join(repoSave, "asterConnectivity.pkl"), "wb") as f:
                    pickle.dump(self._asterConnectivity, f)
        if meshData:
            if self._listPairs is not None:
                with open(os.path.join(repoSave, "listPairs.pkl"), "wb") as f:
                    pickle.dump(self._listPairs, f)
            if self._intePointsList is not None:
                with open(os.path.join(repoSave, "intePointsList.pkl"), "wb") as f:
                    pickle.dump(self._intePointsList, f)
            if self._quadPointsList is not None:
                with open(os.path.join(repoSave, "quadPointsList.pkl"), "wb") as f:
                    pickle.dump(self._quadPointsList, f)
