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

import numpy as np
from libaster import SimpleFieldOnNodesComplex, SimpleFieldOnNodesReal

from ..Objects import PythonBool
from ..Utilities import injector, force_list, medcoupling as medc, ParaMEDMEM as PMM


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

    def toMedCouplingField(self, medmesh):
        """Export the field to a new MEDCoupling field

        Arguments:
            medmesh (*MEDCouplingUMesh*): The medcoupling support mesh.

        Returns:
            field ( MEDCouplingFieldDouble ) : The field medcoupling format.
        """

        if not isinstance(medmesh, (medc.MEDCouplingUMesh, PMM.MEDCouplingUMesh)):
            msg = "toMedCouplingField() argument must be a MEDCouplingUMesh, not '{}'"
            raise TypeError(msg.format(type(medmesh).__name__))

        # Aster values
        values, mask = self.toNumpy()

        # Restrict field based on mask
        restricted_nodes = np.where(np.any(mask, axis=1) == True)[0]
        restricted_values = values[restricted_nodes, :]

        # Med profile
        field_profile = medc.DataArrayInt(restricted_nodes)
        field_profile.setName("NodesProfile")

        # Medcoupling field
        field_values = medc.DataArrayDouble(restricted_values)
        field_values.setInfoOnComponents(self.getComponents())
        medc_node_field = medc.MEDCouplingFieldDouble(medc.ON_NODES, medc.ONE_TIME)
        medc_node_field.setName(self.getPhysicalQuantity())
        medc_node_field.setArray(field_values)
        medc_node_field.setNature(medc.IntensiveMaximum)

        if len(restricted_nodes) == medmesh.getNumberOfNodes():
            medc_node_field.setMesh(medmesh)
        else:
            # Med support mesh for field ( restricted to profile nodes ) without cells
            field_mesh = medc.MEDCouplingUMesh()
            field_mesh.setName("")
            field_mesh.setMeshDimension(medmesh.getMeshDimension())
            field_mesh.setCoords(medmesh.getCoords()[field_profile])
            field_mesh.allocateCells()
            medc_node_field.setMesh(field_mesh)

        medc_node_field.checkConsistencyLight()

        return medc_node_field

    def toMEDFileField1TS(self, medmesh):
        """Export the field to a new MED field

        Arguments:
            medmesh (*MEDFileUMesh*): The medcoupling support mesh.

        Returns:
            field ( MEDFileField1TS ) : The field in med format ( medcoupling ).
        """

        if not isinstance(medmesh, medc.MEDFileUMesh):
            msg = "toMEDFileField1TS() argument must be a MEDFileUMesh, not '{}'"
            raise TypeError(msg.format(type(medmesh).__name__))

        # Aster values
        values, mask = self.toNumpy()

        # Restrict field based on mask
        restricted_nodes = np.where(np.any(mask, axis=1) == True)[0]
        restricted_values = values[restricted_nodes, :]

        # Med profile
        field_profile = medc.DataArrayInt(restricted_nodes)
        field_profile.setName("NodesProfile")

        # Med support mesh for field ( restricted to profile nodes ) without cells
        field_mesh = medc.MEDCouplingUMesh()
        field_mesh.setName("")
        field_mesh.setMeshDimension(medmesh.getMeshDimension())
        field_mesh.setCoords(medmesh.getCoords()[field_profile])
        field_mesh.allocateCells()

        # Medcoupling field
        field_values = medc.DataArrayDouble(restricted_values)
        field_values.setInfoOnComponents(self.getComponents())
        medc_node_field = medc.MEDCouplingFieldDouble(medc.ON_NODES, medc.ONE_TIME)
        medc_node_field.setMesh(field_mesh)
        medc_node_field.setName(self.getPhysicalQuantity())
        medc_node_field.setArray(field_values)
        medc_node_field.checkConsistencyLight()

        # Med field with profile
        medfield = medc.MEDFileField1TS()
        medfield.setFieldProfile(medc_node_field, medmesh, 1, field_profile)

        return medfield


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
