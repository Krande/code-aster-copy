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

# person_in_charge: francesco.bettonte@edf.fr
"""
:py:class:`SimpleFieldOnCellsReal`
Simple Fields defined on cells of elements
********************************************************************
"""

import numpy as np
from libaster import SimpleFieldOnCellsReal

from ..Utilities import injector, medcoupling as medc, ParaMEDMEM as PMM
from ..Objects.Serialization import InternalStateBuilder


class SFoCStateBuilder(InternalStateBuilder):
    """Class that returns the internal state of a *SimpleFieldOnCells*."""

    def restore(self, field):
        """Restore the *DataStructure* content from the previously saved internal
        state.

        Arguments:
            field (*DataStructure*): The *DataStructure* object to be restored.
        """
        super().restore(field)
        field.build()


class CFValue:
    """Represents the values of a component of a field on cells."""

    def __init__(self, values, mask):
        self._values = values
        self._mask = mask

    def size(self):
        """Return the size of the vector."""
        return len(self._values)

    def __repr__(self):
        lines = repr(self._values).splitlines()
        if len(lines) > 6:
            lines = lines[:3] + ["..."] + lines[-3:]
        return "\n".join(lines)

    def _check_consistency(self, other):
        """Check consistency between the both masks. A float is returned as is."""
        if isinstance(other, CFValue):
            if not np.all(self._mask == other._mask):
                raise IndexError("inconsistent description")
            other = other._values
        return other

    def apply(self, func):
        """Apply a function of the values."""
        if not isinstance(func, np.ufunc):
            func = np.vectorize(func)
        return CFValue(func(self._values), self._mask.copy())

    def __add__(self, other):
        other = self._check_consistency(other)
        return CFValue(self._values + other, self._mask.copy())

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + (-other)

    def __neg__(self):
        return CFValue(-self._values, self._mask.copy())

    def __mul__(self, other):
        other = self._check_consistency(other)
        return CFValue(self._values * other, self._mask.copy())

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        other = self._check_consistency(other)
        if np.any(other == 0.0):
            raise ZeroDivisionError()
        return CFValue(self._values / other, self._mask.copy())

    def __rtruediv__(self, other):
        if np.any(self._values == 0.0):
            raise ZeroDivisionError()
        return CFValue(1.0 / self._values * other, self._mask.copy())

    def __abs__(self):
        return CFValue(np.abs(self._values), self._mask.copy())

    def __pow__(self, other):
        other = self._check_consistency(other)
        return CFValue(np.power(self._values, other), self._mask.copy())

    def min(self):
        return self._values.min()

    def max(self):
        return self._values.max()

    def sum(self):
        return self._values.sum()

    def mean(self):
        return self._values.mean()


@injector(SimpleFieldOnCellsReal)
class ExtendedSimpleFieldOnCellsReal:
    internalStateBuilder = SFoCStateBuilder

    @property
    def _cache(self):
        if self._ptr_cache is None:
            self._ptr_cache = dict.fromkeys("vmc")
        if self._ptr_cache["c"] is None:
            self._ptr_cache["c"] = self.getComponents()
        if self._ptr_cache["v"] is None:
            self._ptr_cache["v"], self._ptr_cache["m"] = self.toNumpy()
        return self._ptr_cache

    def __getattr__(self, component):
        idx = self._cache["c"].index(component)
        return CFValue(self._cache["v"][:, idx].copy(), self._cache["m"][:, idx].copy())

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
        values, mask = self._cache["v"], self._cache["m"]
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

    def transfert(self, mesh, cmps=[]):
        """Tranfert the field to an other mesh. One of the mesh has to be a restriction
        to the other one.

        Arguments:
            mesh (Mesh) : mesh to use for transfert.
            cmps [list[str]]: filter on list of components. If empty, all components are used

        Returns:
            SimpleFieldOnCellsReal: field transfered to new mesh.
        """

        if len(cmps) == 0:
            cmps = self.getComponents()
        else:
            cmps_red = []
            all_cmps = self.getComponents()
            for cmp in cmps:
                if cmp in all_cmps:
                    cmps_red.append(cmp)
            cmps = cmps_red

        sfield = SimpleFieldOnCellsReal(
            mesh, self.getLocalization(), self.getPhysicalQuantity(), cmps
        )

        rest2orig = mesh.getRestrictedToOriginalCellsIds()

        if len(rest2orig) > 0:
            orig2rest = mesh.getOriginalToRestrictedCellsIds()

            values, descr = self.getValuesWithDescription(cmps, rest2orig)

            sfield.setValues(
                [orig2rest[cell] for cell in descr[0]], descr[1], descr[2], descr[3], values
            )
        else:
            rest2orig = self.getMesh().getRestrictedToOriginalCellsIds()

            values, descr = self.getValuesWithDescription(cmps)

            sfield.setValues(
                [rest2orig[cell] for cell in descr[0]], descr[1], descr[2], descr[3], values
            )

        return sfield

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

        values, mask = self._cache["v"], self._cache["m"]

        # Restrict field based on mask
        restricted_cells = np.where(np.any(mask, axis=1))[0]
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

    def fromMEDCouplingField(self, mc_field):
        """Import the field to a new MEDCoupling field. Set values in place.

            It assumes that the DataArray contains the list of components and
            its names is the physical quantity.

        Arguments:
            field (*MEDCouplingFieldDouble*): The medcoupling field.
        """

        if not isinstance(mc_field, (medc.MEDCouplingFieldDouble, PMM.MEDCouplingFieldDouble)):
            msg = "fromMEDCouplingField() argument must be a MEDCouplingFieldDouble, not '{}'"
            raise TypeError(msg.format(type(mc_field).__name__))

        if mc_field.getTypeOfField() != medc.ON_CELLS:
            raise RuntimeError("Field is not defined on cells.")

        mesh = mc_field.getMesh()
        if mesh.getNumberOfCells() != self.getMesh().getNumberOfCells():
            raise RuntimeError("Meshes have an incompatible number of cells.")

        # mc values
        arr = mc_field.getArray()

        cmps = arr.getInfoOnComponents()
        phys = arr.getName()

        self.allocate("ELEM", phys, cmps, 1, 1)

        self.setValues(arr.getValues())
