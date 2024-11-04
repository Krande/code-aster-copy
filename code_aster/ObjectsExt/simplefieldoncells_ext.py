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
:py:class:`SimpleFieldOnCellsReal`
Simple Fields defined on cells of elements
********************************************************************
"""

import numpy as np
from libaster import SimpleFieldOnCellsReal

from ..Utilities import injector, medcoupling as medc, ParaMEDMEM as PMM


@injector(SimpleFieldOnCellsReal)
class ExtendedSimpleFieldOnCellsReal:
    def getValues(self, copy=False):
        """
        Returns two numpy arrays containing the field values on specific cells.

        The first array contains the field values while the second one is a mask
        which is `True` if the corresponding value exists, `False` otherwise.
        Array shape is ( number_of_cells_with_components, number_of_components ).

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

    def getValuesOnCell(self, idcell):
        """
        Returns two numpy arrays containing the field values on specific cells.

        The first array contains the field values while the second one is a mask
        which is `True` if the corresponding value exists, `False` otherwise.
        Each array contains `number_of_cells * number_of_components` values.
        Array shape is ( number_of_cells_with_components, number_of_components ).

        Args:
            idcell (int): The cell id.

        Returns:
            ndarray (float): Field values.
            ndarray (bool): Mask for the field values.
        """

        istart = sum(
            [
                self.getNumberOfPointsOfCell(i) * self.getNumberOfSubPointsOfCell(i)
                for i in range(idcell)
            ]
        )

        iend = istart + self.getNumberOfPointsOfCell(idcell) * self.getNumberOfSubPointsOfCell(
            idcell
        )

        values, mask = self.getValues()

        return values[istart:iend], mask[istart:iend]

    def toMEDCouplingField(self, medmesh):
        """Export the field to a new MEDCoupling field

        Arguments:
            medmesh (*MEDCouplingUMesh*): The medcoupling support mesh.

        Returns:
            field ( MEDCouplingFieldDouble ) : The field medcoupling format.
        """

        if not isinstance(medmesh, (medc.MEDCouplingUMesh, PMM.MEDCouplingUMesh)):
            msg = "toMEDCouplingField() argument must be a MEDCouplingUMesh, not '{}'"
            raise TypeError(msg.format(type(medmesh).__name__))

        # Aster values
        values, mask = self.toNumpy()

        # Restrict field based on mask
        restricted_cells = np.where(np.any(mask, axis=1) == True)[0]
        restricted_values = values[restricted_cells, :]

        # Medcoupling field
        field_values = medc.DataArrayDouble(restricted_values)
        field_values.setInfoOnComponents(self.getComponents())
        field_values.setName(self.getPhysicalQuantity())
        medc_cell_field = medc.MEDCouplingFieldDouble(medc.ON_CELLS, medc.ONE_TIME)
        medc_cell_field.setName(self.getPhysicalQuantity())
        medc_cell_field.setArray(field_values)
        medc_cell_field.setNature(medc.IntensiveConservation)

        if len(restricted_cells) == medmesh.getNumberOfCells():
            medc_cell_field.setMesh(medmesh)
        else:
            raise NotImplementedError()

        medc_cell_field.checkConsistencyLight()

        return medc_cell_field
