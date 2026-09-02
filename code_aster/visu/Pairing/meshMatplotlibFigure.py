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
from .convexPointSet import ConvexPointSet
import numpy as np

try:
    import matplotlib

    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

DEFAULT_FIGURE_SIZE = (11, 9)
DEFAULT_TICKS_SIZE = 12
## -----------------------------------------------------------
#   OPTION AVAILABLE FOR THE VISUALISATION
## -----------------------------------------------------------
OPTION_MESH_VISU = {"domain", "interface", "selectSlvCell"}
SUBOPTION_MESH_VISU = {"all", "givenPair", "givenSlvIndex"}
OPTION_PAIRING_VISU = {"meshOnly", "pairs", "intePoints", "quadPoints"}
OPTION_PAIRING_VISU_DETAILS = OPTION_PAIRING_VISU - {"meshOnly"}

DIM_AVAILABLE = {2, 3}

INDEX_PLANE_PROJECTED = {"X", "Y", "Z", None}

VALID_COMBINATION = {
    2: {
        "domain": {
            "all": OPTION_PAIRING_VISU,
            "givenPair": OPTION_PAIRING_VISU_DETAILS,
            "givenSlvIndex": OPTION_PAIRING_VISU_DETAILS,
        },
        "interface": {
            "all": OPTION_PAIRING_VISU,
            "givenPair": OPTION_PAIRING_VISU_DETAILS,
            "givenSlvIndex": OPTION_PAIRING_VISU_DETAILS,
        },
        "selectSlvCell": {"givenSlvIndex": OPTION_PAIRING_VISU_DETAILS},
    },
    3: {
        "domain": {"all": {}, "givenPair": {}, "givenSlvIndex": {}},
        "interface": {
            "all": OPTION_PAIRING_VISU,
            "givenPair": OPTION_PAIRING_VISU_DETAILS,
            "givenSlvIndex": OPTION_PAIRING_VISU_DETAILS,
        },
        "selectSlvCell": {"givenSlvIndex": OPTION_PAIRING_VISU_DETAILS},
    },
}


## -----------------------------------------------------------
#   CLASS MESH MATPLOTLIB FIGURE
## -----------------------------------------------------------
class meshMatplotlibFigure:

    def __init__(
        self,
        pairingAnalysisInstance,
        dimMatPlot,
        optionMesh,
        suboptionMesh,
        optionPair,
        addNodeLabel=False,
        addMeshNodes=False,
        addLegend=True,
        index=None,
        indexPlaneProjected=None,
    ):
        r"""Constructor

        Args:
            pairingAnalysisInstance (:class:`PairingObject`):
            dimMatPlot (:class:`int`): space dimension
            optionMesh (:class:`str`): should be in OPTION_MESH_VISU
            suboptionMesh (:class:`str`): should be in SUBOPTION_MESH_VISU
            optionPair (:class:`str`): should be in OPTION_PAIRING_VISU
            addNodeLabel (:class:`bool`): if True, then add node labels (numbering)
            addMeshNodes (:class:`bool`): if True, then add mesh node (bullet for nodes)
            addLegend (:class:`bool`): if True, then plot legend
            index (:class:`int`), optional: index of a selected cell
            indexPlaneProjected (:class:`str`), optional: should be un INDEX_PLANE_PROJECTED
        """
        self._pairingAnalysis = pairingAnalysisInstance
        # - Option for visualisation
        self._optionMesh = optionMesh
        self._suboptionMesh = suboptionMesh
        self._optionPair = optionPair
        # - Details for the plot
        self._dimMatPlot = dimMatPlot
        self._addNodeLabel = addNodeLabel
        self._addMeshNodes = addMeshNodes
        self._addLegend = addLegend
        self._index = index
        self._indexPlaneProjected = indexPlaneProjected
        # - Initialisation of some variables
        self._dim = None
        self._codim = None
        # - Precomputation
        self.checkConsistency()
        self.computeDimCodim()
        self.setIndicesProjected()
        self.validityOptions()

    def checkConsistency(self):
        r"""Consistency check of the options"""
        if self._optionMesh not in OPTION_MESH_VISU:
            raise ValueError(f"Key {self._optionMesh} not in OPTION_MESH_VISU definition")
        if self._suboptionMesh not in SUBOPTION_MESH_VISU:
            raise ValueError(f"Key {self._suboptionMesh} not in SUBOPTION_MESH_VISU definition")
        if self._optionPair not in OPTION_PAIRING_VISU:
            raise ValueError(f"Key {self._optionPair} not in OPTION_PAIRING_VISU definition")
        if self._indexPlaneProjected not in INDEX_PLANE_PROJECTED:
            raise ValueError(
                f"Key : indexPlaneProjected '{self._indexPlaneProjected}' not in INDEX_PLANE_PROJECTED"
            )

    def validityOptions(self):
        # - Check whether the method is valid for the given dimension
        if self._suboptionMesh not in VALID_COMBINATION[self._dim][self._optionMesh]:
            raise ValueError(
                f"ERROR : Option '{self._suboptionMesh}' invalid for the method '{self._optionMesh}' and dimension '{self._dim}'."
            )
        # - Check if the option is valid for the given method and dimension
        if self._suboptionMesh not in VALID_COMBINATION[self._dim][self._optionMesh]:
            raise ValueError(
                f"ERROR : Option '{self._suboptionMesh}' invalid for the method '{self._optionMesh}' and dimension '{self._dim}'."
            )
        # - Check if the sub-option is valid for the given option
        if (
            self._optionPair
            not in VALID_COMBINATION[self._dim][self._optionMesh][self._suboptionMesh]
        ):
            raise ValueError(
                f"ERROR : Suboption '{self._optionPair}' invalid for option '{self._suboptionMesh}' and method '{self._optionMesh}'."
            )

    def computeDimCodim(self):
        r"""Method to computed codimension"""
        if self._optionMesh == "domain":
            self._dim = self._dimMatPlot
            self._codim = 0
        elif self._optionMesh == "interface":
            self._dim = self._dimMatPlot
            self._codim = 1
        elif self._optionMesh == "selectSlvCell":
            self._dim = self._dimMatPlot
            self._codim = 1
        else:
            raise ValueError(f"Key {self._optionMesh} not in OPTION_MESH_VISU definition")

    def setIndicesProjected(self):
        r"""Given the option provided, defines the coordinates to be selected
        for a 3D plot projected onto a plane."""
        if self._indexPlaneProjected == "X":
            self._indexPPx = 1
            self._indexPPy = 2
            self._indexPPz = 0  # - coordinates removed
        elif self._indexPlaneProjected == "Y":
            self._indexPPx = 0
            self._indexPPy = 2
            self._indexPPz = 1  # - coordinates removed
        elif self._indexPlaneProjected == "Z":
            self._indexPPx = 0
            self._indexPPy = 1
            self._indexPPz = 2  # - coordinates removed
        elif self._indexPlaneProjected is None:
            self._indexPPx = None
            self._indexPPy = None
            self._indexPPz = None  # - coordinates removed
        else:
            raise ValueError(
                f"Error : indexPlaneProjected '{self._indexPlaneProjected}' is not valid."
            )

    def generateAxis(self):
        r"""Axis generation for plot"""
        if self._dimMatPlot == 2:
            fig, ax = plt.subplots(figsize=DEFAULT_FIGURE_SIZE)
            self._fig = fig
            self._ax = ax
            plt.tick_params(axis="both", labelsize=DEFAULT_TICKS_SIZE)
        elif self._dimMatPlot == 3 and self._indexPlaneProjected is None:
            plt.ion()
            self._fig = plt.figure(figsize=DEFAULT_FIGURE_SIZE)
            self._ax = plt.axes(projection="3d")
            plt.tick_params(axis="both", labelsize=DEFAULT_TICKS_SIZE)
        elif self._dimMatPlot == 3 and self._indexPlaneProjected is not None:
            fig, ax = plt.subplots(figsize=DEFAULT_FIGURE_SIZE)
            plt.tick_params(axis="both", labelsize=DEFAULT_TICKS_SIZE)
            self._fig = fig
            self._ax = ax
        else:
            raise NameError("generateAxis: Dimension of plot not available")

    def addNodes(self, nodesCoords, color, s):
        r"""Add nodes in a figure

        Args:
            nodesCoords (:class:`numpy.ndarray`): Array of the nodes coordinates
            color (:class:`str`): Color of the node coordinates
            s (:class:`float`): Opacity parameter

        """
        if self._dimMatPlot == 2:
            self._ax.scatter(nodesCoords[:, 0], nodesCoords[:, 1], color=color, s=s)
        elif self._dimMatPlot == 3 and self._indexPlaneProjected is None:
            self._ax.scatter(
                nodesCoords[:, 0], nodesCoords[:, 1], nodesCoords[:, 2], color=color, s=s
            )
        elif self._dimMatPlot == 3 and self._indexPlaneProjected is not None:
            self._ax.scatter(
                nodesCoords[:, self._indexPPx], nodesCoords[:, self._indexPPy], color=color, s=s
            )
        else:
            raise NameError("addNodes: Dimension of plot not available")

    def addNodeLabels(self, nodeCoords, nodeIndices, color):
        r"""Add node labels to the plot

        Args:
            nodeCoords (:class:`numpy.ndarray`): Array of nodes coordinates
            color (:class:`str`): Color of the edges

        """
        factor = 1.01
        for nodeIndex in range(nodeCoords.shape[0]):
            if self._dimMatPlot == 2:
                self._ax.text(
                    factor * nodeCoords[nodeIndex, 0],
                    factor * nodeCoords[nodeIndex, 1],
                    "%s" % (nodeIndices[nodeIndex]),
                    size=15,
                    zorder=2,
                    color=color,
                )
            elif self._dimMatPlot == 3 and self._indexPlaneProjected is None:
                self._ax.text(
                    factor * nodeCoords[nodeIndex, 0],
                    factor * nodeCoords[nodeIndex, 1],
                    factor * nodeCoords[nodeIndex, 2],
                    "%s" % (nodeIndices[nodeIndex]),
                    size=15,
                    zorder=2,
                    color=color,
                )
            elif self._dimMatPlot == 3 and self._indexPlaneProjected is not None:
                self._ax.text(
                    factor * nodeCoords[nodeIndex, self._indexPPx],
                    factor * nodeCoords[nodeIndex, self._indexPPy],
                    "%s" % (nodeIndices[nodeIndex]),
                    size=15,
                    zorder=2,
                    color=color,
                )
            else:
                raise NameError("addNodeLabels: Dimension of plot not available")

    def addEdges(self, cellIndices, dim, codim, color, linestyle, alpha, labelStr):
        r"""Add edges associate to a cell

        Args:
            cellIndices (:class:`list`): List of the cell indices to plot
            dim (:class:`int`): Dimension of the space
            codim (:class:`int`): Codimension for the convex set to plot
            color (:class:`str`): Color of the edges
            linestyle (:class:`str`): Linestyle of the edges
            alpha (:class:`float`): Opacity parameter
            labelStr (:class:`float`): Label for this cell

        """
        for k, index in enumerate(cellIndices):
            nodes, _ = self._pairingAnalysis.getNodesCoordsFromCellIndices([index])
            if self._indexPlaneProjected is not None:
                indicesProj = [self._indexPPx, self._indexPPy]
                nodesNew = np.zeros(np.shape(nodes))
                nodesExtr = nodes[:, indicesProj]
                nodesNew[:, [0, 1]] = nodesExtr
                cell = ConvexPointSet(dim - codim, nodesNew, 0)
            else:
                cell = ConvexPointSet(dim - codim, nodes, codim)
            if self._addLegend and k == 0:
                label = labelStr
            else:
                label = None

            if self._indexPlaneProjected is not None:
                cell.plotEdges(self._ax, color, linestyle, alpha, self._dim - 1, label)
            else:
                cell.plotEdges(self._ax, color, linestyle, alpha, self._dim, label)

    def setCellsIndices(self):
        r"""Set indices of the slave and master cells"""
        if self._optionMesh == "domain":
            self._SlvIndices = self._pairingAnalysis._indicesSlaveDomain
            self._MasIndices = self._pairingAnalysis._indicesMasterDomain
            # codim = 0
        elif self._optionMesh == "interface":
            self._SlvIndices = self._pairingAnalysis._indicesSlaveInterface
            self._MasIndices = self._pairingAnalysis._indicesMasterInterface
            # codim = 1
        elif self._optionMesh == "selectSlvCell":
            self._SlvIndices = self._pairingAnalysis._indicesSlaveInterface
            self._MasIndices = self._pairingAnalysis._indicesMasterInterface
            # codim = 1
        else:
            raise ValueError(f"Key {self._optionMesh} not in OPTION_MESH_VISU definition")

    def setLegend(self):
        r"""Add legend to the matplotlib figure"""
        if self._addLegend:
            if self._optionPair:
                ncolLeg = 4
            else:
                ncolLeg = 2
            self._ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.1), ncol=ncolLeg)

    def plotMeshStructure(self, plotParams, s):
        r"""Plot the mesh structure (nodes and edges

        Args:
            plotParams (:class:`dict`): dictionnary of parameters for plot (colors)
            s (:class:`float`): Opacity parameter
        """
        # - Slave and Master get indices
        nodesSlvcoords, nodesSlvIndices = self._pairingAnalysis.getNodesCoordsFromCellIndices(
            self._SlvIndices
        )
        nodesMascoords, nodesMasIndices = self._pairingAnalysis.getNodesCoordsFromCellIndices(
            self._MasIndices
        )
        # - Add mesh nodes and/or labels if needed
        if self._addMeshNodes:
            # - Slave
            self.addNodes(nodesSlvcoords, plotParams["Slv"]["color"], s)
            if self._addNodeLabel:
                for nodeIndex in range(nodesSlvcoords.shape[0]):
                    self.addNodeLabels(nodesSlvcoords, nodesSlvIndices, plotParams["Slv"]["color"])
            # - Master
            self.addNodes(nodesMascoords, plotParams["Mas"]["color"], s)
            if self._addNodeLabel:
                for nodeIndex in range(nodesMascoords.shape[0]):
                    self.addNodeLabels(nodesMascoords, nodesMasIndices, plotParams["Mas"]["color"])
        # - Slave and Master edges plot
        self.addEdges(
            self._SlvIndices,
            self._dim,
            self._codim,
            plotParams["Slv"]["color"],
            plotParams["Slv"]["linestyle"],
            plotParams["Slv"]["alpha"],
            "slave",
        )
        self.addEdges(
            self._MasIndices,
            self._dim,
            self._codim,
            plotParams["Mas"]["color"],
            plotParams["Mas"]["linestyle"],
            plotParams["Mas"]["alpha"],
            "master",
        )

    def plotPair(
        self, pairSlvInd, pairMasInd, pairIndex, plotParams, intePts=False, inteQuad=False
    ):
        r"""Plot a pair

        Args:
            pairSlvInd (:class:`int`): index of the slave cell
            pairMasInd (:class:`int`): index of the master cell
            pairIndex (:class:`int`): pair index
            plotParams (:class:`dict`): plot parameters
            intePts (:class:`bool`):
            inteQuad (:class:`bool`):
        """
        if intePts == False:
            if self._indexPlaneProjected is None:
                # - Slave pair
                nodes, _ = self._pairingAnalysis.getNodesCoordsFromCellIndices([pairSlvInd])
                cell = ConvexPointSet(self._dim - self._codim, nodes, self._codim)
                cell.plotEdges(
                    self._ax,
                    plotParams["SlvInt"]["color"],
                    plotParams["SlvInt"]["linestyle"],
                    plotParams["SlvInt"]["alpha"],
                    self._dim,
                    "cell=" + str(pairSlvInd),
                )
                # - Master pair
                nodes, _ = self._pairingAnalysis.getNodesCoordsFromCellIndices([pairMasInd])
                cell = ConvexPointSet(self._dim - self._codim, nodes, self._codim)
                cell.plotEdges(
                    self._ax,
                    plotParams["MasInt"]["color"],
                    plotParams["MasInt"]["linestyle"],
                    plotParams["MasInt"]["alpha"],
                    self._dim,
                    "cell=" + str(pairMasInd),
                )
            else:
                # - Slave pair
                nodes, _ = self._pairingAnalysis.getNodesCoordsFromCellIndices([pairSlvInd])
                indicesProj = [self._indexPPx, self._indexPPy]
                nodesNew = np.zeros(np.shape(nodes))
                nodesExtr = nodes[:, indicesProj]
                nodesNew[:, [0, 1]] = nodesExtr
                cell = ConvexPointSet(self._dim - self._codim, nodesNew, 0)
                # cell = ConvexPointSet(self.dim - self.codim, nodes, self.codim)
                cell.plotEdges(
                    self._ax,
                    plotParams["SlvInt"]["color"],
                    plotParams["SlvInt"]["linestyle"],
                    plotParams["SlvInt"]["alpha"],
                    self._dim - 1,
                    "cell=" + str(pairSlvInd),
                )
                # - Master pair
                nodes, _ = self._pairingAnalysis.getNodesCoordsFromCellIndices([pairMasInd])
                indicesProj = [self._indexPPx, self._indexPPy]
                nodesNew = np.zeros(np.shape(nodes))
                nodesExtr = nodes[:, indicesProj]
                nodesNew[:, [0, 1]] = nodesExtr
                cell = ConvexPointSet(self._dim - self._codim, nodesNew, 0)
                # cell = ConvexPointSet(self.dim - self.codim, nodes, self.codim)
                cell.plotEdges(
                    self._ax,
                    plotParams["MasInt"]["color"],
                    plotParams["MasInt"]["linestyle"],
                    plotParams["MasInt"]["alpha"],
                    self._dim - 1,
                    "cell=" + str(pairMasInd),
                )
            # raise ValueError("not implemented")
        else:
            if inteQuad == False:
                raise ValueError("not implemented")
            else:
                raise ValueError("not implemented")

    def plotIntePts(self, pairIndex, s):
        r"""Plot the intersection points for a given pair.

        Args:
            pairIndex (:class:`int`): index of the pair in the list of pairs
            s (:class:`float`): Opacity parameter
        """
        if self._indexPlaneProjected is None:
            # - Intersection points
            inteConvexSet = [
                list(tu) for tu in self._pairingAnalysis._listIntersectionPts[pairIndex]
            ]
            inteConvexSetNP = np.array(inteConvexSet)
            cell = ConvexPointSet(self._dim - self._codim, inteConvexSet, self._codim)
            cell.plotEdges(self._ax, "black", "-", 1.0, self._dim, None)
        else:
            # - Intersection points
            inteConvexSet = [
                list(tu) for tu in self._pairingAnalysis._listIntersectionPts[pairIndex]
            ]
            inteConvexSetNP = np.array(inteConvexSet)
            indicesProj = [self._indexPPx, self._indexPPy]
            inteConvexSetNPExtr = inteConvexSetNP[:, indicesProj]
            inteConvexSetNPNew = np.zeros(np.shape(inteConvexSetNP))
            inteConvexSetNPNew[:, [0, 1]] = inteConvexSetNPExtr
            cell = ConvexPointSet(self._dim - self._codim, inteConvexSetNPNew, 0)
            cell.plotEdges(self._ax, "black", "-", 1.0, self._dim - 1, None)
        if self._dimMatPlot == 2:
            self._ax.scatter(inteConvexSetNP[:, 0], inteConvexSetNP[:, 1], color="black", s=s)
        elif self._dimMatPlot == 3 and self._indexPlaneProjected is None:
            self._ax.scatter(
                inteConvexSetNP[:, 0],
                inteConvexSetNP[:, 1],
                inteConvexSetNP[:, 2],
                color="black",
                s=s,
            )
        elif self._dimMatPlot == 3 and self._indexPlaneProjected is not None:
            self._ax.scatter(
                inteConvexSetNP[:, self._indexPPx],
                inteConvexSetNP[:, self._indexPPy],
                color="black",
                s=s,
            )
        else:
            raise ValueError("dimension of the plot is either 2 or 3")

    def plotQuadPts(self, pairIndex, s):
        r"""Plot the quadrature points for a given pair.

        Args:
            pairIndex (:class:`int`): index of the pair in the list of pairs
            s (:class:`float`): Opacity parameter
        """
        if self._indexPlaneProjected is None:
            inteConvexSet = self._pairingAnalysis._listIntersectionPts[pairIndex]
            cell = ConvexPointSet(self._dim - self._codim, inteConvexSet, self._codim)
            cell.computeBary()
            cell.plotEdges(self._ax, "black", "-", 1.0, self._dim, None)
        else:
            inteConvexSet = self._pairingAnalysis._listIntersectionPts[pairIndex]
            indicesProj = [self._indexPPx, self._indexPPy]
            inteConvexSetNew = np.zeros(np.shape(inteConvexSet))
            inteConvexSetExtr = inteConvexSet[:, indicesProj]
            inteConvexSetNew[:, [0, 1]] = inteConvexSetExtr
            cell = ConvexPointSet(self._dim - self._codim, inteConvexSetNew, 0)
            cell.computeBary()
            cell.plotEdges(self._ax, "black", "-", 1.0, self._dim - 1, None)
        for k in range(len(inteConvexSet)):
            if self._dimMatPlot == 2:
                plt.plot(
                    [inteConvexSet[k][0], cell.computeBary()[0]],
                    [inteConvexSet[k][1], cell.computeBary()[1]],
                    color="black",
                    marker=None,
                    linestyle="dashed",
                )
            elif self._dimMatPlot == 3 and self._indexPlaneProjected is None:
                plt.plot(
                    [inteConvexSet[k][0], cell.computeBary()[0]],
                    [inteConvexSet[k][1], cell.computeBary()[1]],
                    [inteConvexSet[k][2], cell.computeBary()[2]],
                    color="black",
                    marker=None,
                    linestyle="dashed",
                )
            elif self._dimMatPlot == 3 and self._indexPlaneProjected is not None:
                plt.plot(
                    [inteConvexSet[k][self._indexPPx], cell.computeBary()[self._indexPPx]],
                    [inteConvexSet[k][self._indexPPy], cell.computeBary()[self._indexPPy]],
                    color="black",
                    marker=None,
                    linestyle="dashed",
                )
            else:
                raise ValueError("dimension of the plot is either 2 or 3")
        # - Quadrature points
        quadPointspairIndex = self._pairingAnalysis._listQuadraturePts[pairIndex]
        quadPointspairIndexNP = np.array(quadPointspairIndex)
        if self._dimMatPlot == 2:
            self._ax.scatter(
                quadPointspairIndexNP[:, 0],
                quadPointspairIndexNP[:, 1],
                color="black",
                s=s,
                marker="o",
            )
        elif self._dimMatPlot == 3 and self._indexPlaneProjected is None:
            self._ax.scatter(
                quadPointspairIndexNP[:, 0],
                quadPointspairIndexNP[:, 1],
                quadPointspairIndexNP[:, 2],
                color="black",
                s=s,
                marker="o",
            )
        elif self._dimMatPlot == 3 and self._indexPlaneProjected is not None:
            self._ax.scatter(
                quadPointspairIndexNP[:, self._indexPPx],
                quadPointspairIndexNP[:, self._indexPPy],
                color="black",
                s=s,
                marker="o",
            )
        else:
            raise ValueError("dimension of the plot is either 2 or 3")

    def plot(self, s=50):
        r"""Method for generating a plot given pairing coordinates"""
        if HAS_MATPLOTLIB:
            # - Initialisation of the cells indices to plot
            self.setCellsIndices()
            if self._optionMesh != "selectSlvCell":
                # - Set the plotting parameters
                if self._optionPair == "meshOnly":
                    plotParams = {
                        "Slv": {"color": "blue", "linestyle": "-", "alpha": 1.0},
                        "Mas": {"color": "red", "linestyle": "-", "alpha": 1.0},
                    }
                else:
                    plotParams = {
                        "Slv": {"color": "blue", "linestyle": "-", "alpha": 0.15},
                        "Mas": {"color": "red", "linestyle": "-", "alpha": 0.15},
                        "SlvInt": {"color": "orange", "linestyle": "-", "alpha": 1.0},
                        "MasInt": {"color": "green", "linestyle": "-", "alpha": 1.0},
                    }

                if self._optionPair == "meshOnly":

                    # - Initialisation of the Figure
                    self.generateAxis()
                    # - Plot Mesh
                    self.plotMeshStructure(plotParams, s)
                    self.setLegend()
                    plt.show(block=True)

                # elif self.optionPair == "Pairs":
                else:
                    if self._suboptionMesh == "all":

                        for pairIndex in range(len(self._pairingAnalysis._listPairs)):
                            self.generateAxis()
                            self.plotMeshStructure(plotParams, s)
                            pairSlvInd, pairMasInd = self._pairingAnalysis._listPairs[pairIndex]
                            self.plotPair(
                                pairSlvInd, pairMasInd, pairIndex, plotParams, False, False
                            )
                            if self._optionPair == "intePoints":
                                self.plotIntePts(pairIndex, s)
                            if self._optionPair == "quadPoints":
                                self.plotQuadPts(pairIndex, int(s / 3))
                            self.setLegend()
                            plt.show(block=True)

                    elif self._suboptionMesh == "givenPair":
                        self.generateAxis()
                        self.plotMeshStructure(plotParams, s)
                        pairSlvInd, pairMasInd = self._pairingAnalysis._listPairs[self._index]
                        self.plotPair(pairSlvInd, pairMasInd, self._index, plotParams, False, False)
                        if self._optionPair == "intePoints":
                            self.plotIntePts(self._index, s)
                        if self._optionPair == "quadPoints":
                            self.plotQuadPts(self._index, int(s / 3))
                        self.setLegend()
                        plt.show(block=True)

                    elif self._suboptionMesh == "givenSlvIndex":
                        for pairIndex in range(len(self._pairingAnalysis._listPairs)):
                            pairSlvInd, pairMasInd = self._pairingAnalysis._listPairs[pairIndex]
                            if pairSlvInd == self._index:
                                self.generateAxis()
                                self.plotMeshStructure(plotParams, s)
                                self.plotPair(
                                    pairSlvInd, pairMasInd, pairIndex, plotParams, False, False
                                )
                                if self._optionPair == "intePoints":
                                    self.plotIntePts(pairIndex, s)
                                if self._optionPair == "quadPoints":
                                    self.plotQuadPts(pairIndex, int(s / 3))
                                self.setLegend()
                                plt.show(block=True)
                    else:
                        raise ValueError("not implemented yet")
                # else:
                #     raise ValueError("not implemented yet")
            else:
                plotParams = {
                    "Slv": {"color": "blue", "linestyle": "-", "alpha": 0.15},
                    "Mas": {"color": "red", "linestyle": "-", "alpha": 0.15},
                    "SlvInt": {"color": "blue", "linestyle": "-", "alpha": 1.0},
                    "MasInt": {"color": "red", "linestyle": "-", "alpha": 1.0},
                }
                self.generateAxis()
                for pairIndex in range(len(self._pairingAnalysis._listPairs)):
                    pairSlvInd, pairMasInd = self._pairingAnalysis._listPairs[pairIndex]
                    if pairSlvInd == self._index:
                        # - Slave and Master get indices
                        nodesSlvcoords, _ = self._pairingAnalysis.getNodesCoordsFromCellIndices(
                            [pairSlvInd]
                        )
                        nodesMascoords, _ = self._pairingAnalysis.getNodesCoordsFromCellIndices(
                            [pairMasInd]
                        )
                        if self._addMeshNodes:
                            # - Slave
                            self.addNodes(nodesSlvcoords, plotParams["Slv"]["color"], s)
                            if self._addNodeLabel:
                                for nodeIndex in range(nodesSlvcoords.shape[0]):
                                    self.addNodeLabels(
                                        nodesSlvcoords, [pairSlvInd], plotParams["Slv"]["color"]
                                    )
                            # - Master
                            self.addNodes(nodesMascoords, plotParams["Mas"]["color"], s)
                            if self._addNodeLabel:
                                for nodeIndex in range(nodesMascoords.shape[0]):
                                    self.addNodeLabels(
                                        nodesMascoords, [pairMasInd], plotParams["Mas"]["color"]
                                    )
                        self.plotPair(pairSlvInd, pairMasInd, pairIndex, plotParams, False, False)
                for pairIndex in range(len(self._pairingAnalysis._listPairs)):
                    pairSlvInd, pairMasInd = self._pairingAnalysis._listPairs[pairIndex]
                    if pairSlvInd == self._index:
                        if self._optionPair == "intePoints":
                            self.plotIntePts(pairIndex, s)
                        if self._optionPair == "quadPoints":
                            self.plotQuadPts(pairIndex, int(s / 3))
                self.setLegend()
                plt.show(block=True)
        else:
            raise ValueError(
                f"HAS_MATPLOTLIB value is {HAS_MATPLOTLIB}. Should try to run it in interactive mode."
            )
