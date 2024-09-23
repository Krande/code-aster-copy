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

# person_in_charge: francesco.bettonte@edf.fr
"""
:py:class:`SimpleFieldOnNodesReal`
Simple Fields defined on nodes of elements
********************************************************************
"""

import aster
import numpy as np
from libaster import SimpleFieldOnNodesComplex, SimpleFieldOnNodesReal

from ..Objects import PythonBool
from ..Utilities import injector, force_list, medcoupling as medc


@injector(SimpleFieldOnNodesReal)
class ExtendedSimpleFieldOnNodesReal:
    def restrict(self, cmps=[], groupsOfNodes=[], same_rank=None):
        """Return a new field restricted to the list of components and groups of nodes given

        Arguments:
            cmps[list[str]]: filter on list of components
            If empty, all components are used
            groupsOfNodes[list[str]]: filter on list of groups of nodes (default=" ").
            If empty, the full mesh is used
            same_rank : - None: keep all nodes (default: None)
                        - True: keep the nodes which are owned by the current MPI-rank
                        - False: keep the nodes which are not owned by the current MPI-rank

        Returns:
            SimpleFieldOnNodesReal: field restricted.
        """

        val = {None: PythonBool.NONE, True: PythonBool.TRUE, False: PythonBool.FALSE}

        return self._restrict(force_list(cmps), force_list(groupsOfNodes), val[same_rank])

    def getValues(self, copy=False):
        """
        Returns two numpy arrays containing the field values on specific nodes.
        The first array contains the field values while the second one is a mask
        which is `True` if the corresponding value exists, `False` otherwise.
        Array shape is ( number_of_nodes, number_of_components ).

        Args:
            copy (bool): If True copy the data, default: *False*

        Returns:
            ndarray (float): Field values.
            ndarray (bool): Mask for the field values.
        """
        values, mask = self.toNumpy()
        if copy:
            return values.copy(), mask.copy()

        values.setflags(write=False)
        mask.setflags(write=False)

        return values, mask


@injector(SimpleFieldOnNodesComplex)
class ExtendedSimpleFieldOnNodesComplex:
    def getValues(self, copy=False):
        """
        Returns two numpy arrays with shape ( number_of_cells_with_components, number_of_components )
        The first array contains the field values while the second one is a mask
        which is `True` if the corresponding value exists, `False` otherwise.

        Where the mask is `False` the corresponding value is set to zero.

        Args:
            copy (bool): If True copy the data, default: *False*

        Returns:
            ndarray (float): Field values.
            ndarray (bool): Mask for the field values.
        """
        values, mask = self.toNumpy()
        if copy:
            return values.copy(), mask.copy()

        values.setflags(write=False)
        mask.setflags(write=False)

        return values, mask
