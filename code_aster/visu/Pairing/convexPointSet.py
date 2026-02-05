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


def transformToPlanarCoords(point, ref, vec1, vec2):
    r"""Transforms a 3D point into planar coordinates based on a basis defined by two vectors.
    This function projects the vector from a reference point (ref) to a given point (point)
    onto the two basis vectors (vec1, vec2), and returns the resulting 2D coordinates.

    Args:
        point (:class:`numpy.ndarray`): 3D point to transform (array of 3 elements).
        ref (:class:`numpy.ndarray`): 3D reference point (array of 3 elements).
        vec1 (:class:`numpy.ndarray`): First basis vector (array of 3 elements).
        vec2 (:class:`numpy.ndarray`): Second basis vector (array of 3 elements).

    Returns:
        (:class:`numpy.ndarray`): Array of 2 elements containing the planar coordinates.
    """
    # - Compute the vector between a given point and the ref point (ref)
    pointVector = point - ref
    # - Project on the basis vectors vec1 and vec2
    coord1 = np.dot(pointVector, vec1)
    coord2 = np.dot(pointVector, vec2)

    return np.array([coord1, coord2])


## -----------------------------------------------------------
#   CONVEX SET POINT CLASS
## -----------------------------------------------------------
class ConvexPointSet:
    def __init__(self, dim, points, codim):
        r"""Constructor

        Args:
            dim (:class:`int`): Dimension of the problem.
            points (:class:`list`): List of points to define the convex set
            codim (:class:`int`): Codimension of the convex point set considered

        """
        self._dim = dim
        self._initPoints = np.copy(points)
        assert self._dim in [1, 2, 3]
        if self._dim == 1:
            assert len(self._initPoints) == 2
        self._pointSorted = False
        self._codim = codim

    def computeBary(self):
        r"""Compute the barycenter of the convex set of points"""
        self._bary = np.mean(self._initPoints, axis=0)

    def computeLocalBasis(self):
        r"""Compute the local basis of the convex set of points"""
        if self._dim == 2:
            # - Get initial points (at least 3 in )
            a = self._initPoints[0]
            b = self._initPoints[1]
            c = self._initPoints[2]
            # - Compute two vectors
            vecA = b - a
            vecB = c - a
            if np.allclose(np.cross(vecA, vecB), [0.0, 0.0, 0.0]):
                if self._initPoints.shape[0] == 3:
                    raise NameError("Error: Element is flat!")
                else:  # - BAD AND WILL BREAK ONE DAY
                    vecB = self._initPoints[3] - a
            # - Compute the normal
            normal = np.cross(vecA, vecB)
            norm = normal / np.linalg.norm(normal)
            # - Compute the tangenet vector and normalize them
            vec1 = vecA
            vec2 = np.cross(normal, vec1)
            vec1 = vec1 / np.linalg.norm(vec1)
            vec2 = vec2 / np.linalg.norm(vec2)
            return np.array([vec1, vec2, norm])
        else:
            raise NameError("Use of computeLocalBasis for other dim than 2, not implemented!")

    def sortPoints(self):
        r"""Sort of the points in the convex set of points"""
        if self._dim == 1:
            self._sortedPoints = np.copy(self._initPoints)
            self._pointSorted = True
        elif self._dim == 2:
            if self._codim == 0:
                self.computeBary()
                # - Compute angle of each point according to the barycenter
                angles = np.arctan2(
                    self._initPoints[:, 1] - self._bary[1], self._initPoints[:, 0] - self._bary[0]
                )
                # - Sort points according to this angle
                indicesSorted = np.argsort(angles)
                self._sortedPoints = np.copy(self._initPoints[indicesSorted])
                self._pointSorted = True
            elif self._codim == 1:
                localBasis = self.computeLocalBasis()
                coords = np.zeros((self._initPoints.shape[0], 2))
                for i, point in enumerate(self._initPoints):
                    coordLoc = transformToPlanarCoords(
                        point, self._initPoints[0], localBasis[0], localBasis[1]
                    )
                    coords[i, :] = np.copy(coordLoc)
                self.computeBary()
                baryLoc = transformToPlanarCoords(
                    self._bary, self._initPoints[0], localBasis[0], localBasis[1]
                )
                angles = np.arctan2(coords[:, 1] - baryLoc[1], coords[:, 0] - baryLoc[0])
                # - Sort points according to this angle
                indicesSorted = np.argsort(angles)
                self._sortedPoints = np.copy(self._initPoints[indicesSorted])
                self._pointSorted = True
            else:
                raise NameError("Situation impossible!!")
        else:
            raise NameError("Use of sorting point for convex set in dim >=3 not implemented")

    def computeArea1D(self):
        r"""Compute area of a convex set of points in 1D"""
        return np.linalg.norm(np.subtract(self._initPoints[0, :], self._initPoints[1, :]))

    def computeArea2D(self):
        r"""Compute area of a convex set of points in 2D"""
        n = len(self._sortedPoints)
        area = 0.0
        for i in range(n):
            j = (i + 1) % n  # next point (loop)
            area += (
                self._sortedPoints[i, 0] * self._sortedPoints[j, 1]
                - self._sortedPoints[j, 0] * self._sortedPoints[i, 1]
            )
        return abs(area) / 2.0

    def computeArea(self):
        r"""Compute area of a convex set of points"""
        if not self._pointSorted:
            self.sortPoints()

        if self._dim == 1:
            return self.computeArea1D()
        elif self._dim == 2:
            return self.computeArea2D()
        else:
            raise NameError("computeArea not implemented for convex set in dim >=3")

    def plotEdges(self, axis_, colorGiven, linestyleGiven, alphaGiven, dimSpace, label=None):
        r"""
        Plots the edges of the object as a closed polygon on the given axis.

        This method sorts the points if not already sorted, then plots them as a closed
        polygon (connecting the last point back to the first) on the specified axis.
        The plot is adapted to the dimensionality of the space (2D or 3D).

        Args:
            axis_ (:class:`matplotlib.axes.Axes`):  The axis on which to plot the edges.
            colorGiven  (:class:`str`): The color of the plotted edges.
            linestyleGiven  (:class:`str`):  The linestyle of the plotted edges (e.g., '-', '--', ':').
            alphaGiven  (:class:`str`): The transparency level of the plotted edges (0.0 to 1.0).
            dimSpace  (:class:`str`): The dimensionality of the space (2 for 2D, 3 for 3D).
            label  (:class:`str`), optional: The label for the plot legend. If None, no label is set.
        """
        if not self._pointSorted:
            self.sortPoints()
        initP = np.array([self._sortedPoints[0, :]])
        pointsPlot = np.concatenate((self._sortedPoints, initP), axis=0)
        if dimSpace < 3:
            if label is not None:
                axis_.plot(
                    pointsPlot[:, 0],
                    pointsPlot[:, 1],
                    color=colorGiven,
                    linestyle=linestyleGiven,
                    alpha=alphaGiven,
                    label=label,
                )
            else:
                axis_.plot(
                    pointsPlot[:, 0],
                    pointsPlot[:, 1],
                    color=colorGiven,
                    linestyle=linestyleGiven,
                    alpha=alphaGiven,
                )
        else:
            if label is not None:
                axis_.plot(
                    pointsPlot[:, 0],
                    pointsPlot[:, 1],
                    pointsPlot[:, 2],
                    color=colorGiven,
                    linestyle=linestyleGiven,
                    alpha=alphaGiven,
                    label=label,
                )
            else:
                axis_.plot(
                    pointsPlot[:, 0],
                    pointsPlot[:, 1],
                    pointsPlot[:, 2],
                    color=colorGiven,
                    linestyle=linestyleGiven,
                    alpha=alphaGiven,
                )
