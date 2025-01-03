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

from ..Utilities import ParaMEDMEM as PMM
from ..Utilities import force_list, injector, logger, no_new_attributes
from ..Utilities import medcoupling as medc
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


class ComponentOnCells:
    r"""Represents the values of a component of a field on cells.

    The input arrays will not be modified. The caller should not change them
    during the existence of the *ComponentOnCells* object.

    Args:
        values (ndarray[float]): Values of the component,
            dimension \Sum_i(nbcells_i * nbval_i).
        descr (*Description*): Description of the component.
    """

    class Description:
        """Internal description of a component.

        Args:
            indexes (ndarray[int]): Indexes in *values* by cells and points on cells.
            discr (ndarray[int]): Local dicretization: number of points and subpoints
                by cell.
            mask (ndarray[bool]): Indicator that tells a value exist, same dimension.
        """

        _idx = _mask = _discr = _cells = _nbcells = _nbval = None
        __setattr__ = no_new_attributes(object.__setattr__)

        def __init__(self, indexes, discr, mask):
            self._idx = indexes
            self._discr = discr
            self._mask = mask
            self._cells = None
            # number of cells in this object before restriction
            self._nbcells = indexes.shape[0]
            # number of values stored for the component before any restriction
            self._nbval = None  # = indexes.max() + 1 may be greater, assigned with values.

        def copy(self):
            """Copy of the current description.

            Returns:
                *Decription*: copy of the current description.
            """
            new = __class__(self._idx, self._discr, self._mask)
            new._cells = self._cells
            new._nbcells = self._nbcells
            new._nbval = self._nbval
            return new

        def restrict(self, cells_ids):
            """Restrict the description on given cells.

            Args:
                cells_ids (list[int]): sublist of the existing cells.

            Returns:
                ndarray[int]: Indexes restricted on the cells list.
            """
            if isinstance(cells_ids, int):
                cells_ids = force_list(cells_ids)
            idx = self._idx[cells_ids].ravel()
            idx = idx[idx >= 0]
            self._idx = self._idx[cells_ids]
            self._discr = self._discr[cells_ids]
            self._mask = self._mask[idx]
            self._cells = np.array(cells_ids, dtype=int)
            return idx

        def expand(self, values):
            """Expand the description to the initial shape.

            Args:
                values (ndarray[float]): values to be accordingly expanded.

            Returns:
                *Description*: expanded description.
            """
            idx = self._idx
            idx = idx[idx >= 0]
            new_values = np.zeros(self._nbval, dtype=float)
            new_values[idx] = values
            new_idx = np.ones((self._nbcells, self._idx.shape[1]), dtype=int) * -1
            new_idx[self._cells] = self._idx
            new_discr = np.zeros((self._nbcells, 2), dtype=int)
            new_discr[self._cells] = self._discr
            new_mask = np.zeros(self._nbval, dtype=bool)
            new_mask[idx] = self._mask
            return new_values, ComponentOnCells.Description(new_idx, new_discr, new_mask)

    # WARNING: _idx must be used to access values only a not restricted component.
    # FIXME '_mask' toujours True ?

    _values = _descr = _sign = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, values, description: Description):
        self.init_from(values, description)

    def init_from(self, values, description: Description):
        """Initialize content."""
        self._values = values
        self._descr = description
        self._descr._nbval = len(values)
        self._sign = None

    @property
    def _idx(self):
        return self._descr._idx

    @property
    def _cells(self):
        return self._descr._cells

    def copy(self):
        """Return a copy of the component."""
        return __class__(self._values.copy(), self._descr.copy())

    def getValues(self):
        """Returns the vector values."""
        return self._values

    def getValuesByCells(self):
        """Returns the values by cells, expanded by 0.0 where undefined."""
        values = np.append(self._values, [0.0])
        return np.take(values, self._idx)

    def setValues(self, values):
        """Set the vector values.

        Args:
            values (ndarray[float]): New values.
        """
        assert values.shape == self._values.shape
        self._values = values

    def getCells(self):
        """Return the cells id."""
        if self._cells is None:
            return np.arange(self._descr._nbcells)
        return self._cells

    @property
    def size(self):
        """Return the size of the vector."""
        return len(self._values)

    def __repr__(self):
        lines = repr(self._values).splitlines()
        if len(lines) > 6:
            lines = lines[:3] + ["..."] + lines[-3:]
        return "\n".join(lines)

    @property
    def sign(self):
        """int: Attribute that holds the signature of the component."""
        if self._sign is None:
            self._sign = (self._idx.ravel() * np.arange(self._idx.size)).sum()
        return self._sign

    def _check_consistency(self, other):
        """Check consistency between the both signatures. A float is returned as is."""
        if isinstance(other, ComponentOnCells):
            if self.sign != other.sign:
                logger.warning("lhs: idx.shape: %s, %s", self._idx.shape, self._idx[:5])
                logger.warning("rhs: idx.shape: %s, %s", other._idx.shape, other._idx[:5])
                raise IndexError("inconsistent description")
            other = other._values
        return other

    def apply(self, func):
        """Apply a function of the values."""
        if not isinstance(func, np.ufunc):
            func = np.vectorize(func)
        return ComponentOnCells(func(self._values), self._descr)

    def __add__(self, other):
        other = self._check_consistency(other)
        return ComponentOnCells(self._values + other, self._descr)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + (-other)

    def __neg__(self):
        return ComponentOnCells(-self._values, self._descr)

    def __mul__(self, other):
        other = self._check_consistency(other)
        return ComponentOnCells(self._values * other, self._descr)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        other = self._check_consistency(other)
        if np.any(other == 0.0):
            raise ZeroDivisionError()
        return ComponentOnCells(self._values / other, self._descr)

    def __rtruediv__(self, other):
        if np.any(self._values == 0.0):
            raise ZeroDivisionError()
        return ComponentOnCells(1.0 / self._values * other, self._descr)

    def __abs__(self):
        return ComponentOnCells(np.abs(self._values), self._descr)

    def __pow__(self, other):
        other = self._check_consistency(other)
        return ComponentOnCells(np.power(self._values, other), self._descr)

    def min(self):
        return self._values.min()

    def max(self):
        return self._values.max()

    def sum(self):
        return self._values.sum()

    def mean(self):
        return self._values.mean()

    def restrict(self, cells_ids):
        """Restrict the component on given cells."""
        if self._cells is not None:
            exp = self.expand()
            self.init_from(exp._values, exp._descr)
        idx = self._descr.restrict(cells_ids)
        self._values = self._values[idx]

    def expand(self):
        """Expand the component by adding 0. where it not defined.

        Returns:
            *ComponentOnCells*: expanded component.
        """
        if self._cells is None:
            return self.copy()
        return ComponentOnCells(*self._descr.expand(self._values))

    def on_support_of(self, other, strict=False):
        """Move the component onto the same support of another.

        NB: The both components must be defined on the same cells.

        Args:
            other (*ComponentOnCells*): Use the support of this component.
            strict (bool): If *True* it raises an *IndexError* exception if
                the component is not defined on the support of `other`.
                By default, it is *True* and the values are expanded with 0.0.

        Returns:
            *ComponentOnCells*: A new component.
        """
        if other._cells is not None and self._cells is None:
            same = self.copy()
            same.restrict(other._cells)
            return same.on_support_of(other)
        if self._cells is not None and other._cells is None:
            other = other.copy()
            other.restrict(self._cells)
            return self.on_support_of(other)
        if self._cells is not None:
            if np.any(self._cells != other._cells):
                raise IndexError("components must be defined on the same cells")
            homo = self.expand().on_support_of(other.expand())
            if self._cells is not None:
                homo.restrict(self._cells)
            return homo
        new = other.copy()
        keep = self._idx[other._idx >= 0]
        if strict and np.any(keep < 0):
            raise IndexError("no value on all the support")
        prol = np.append(self._values, [0.0])
        new._values = np.take(prol, keep)
        new._descr._idx = other._idx
        return new


@injector(SimpleFieldOnCellsReal)
class ExtendedSimpleFieldOnCellsReal:
    internalStateBuilder = SFoCStateBuilder

    @property
    def _cache(self):
        if self._ptr_cache is None:
            self._ptr_cache = dict.fromkeys(["readonly", "val", "msk", "cmp", "idx", "nbpt"])
            self._ptr_cache["readonly"] = None
        if self._ptr_cache["cmp"] is None:
            self._ptr_cache["cmp"] = self.getComponents()
        if self._ptr_cache["val"] is None:
            self._ptr_cache["val"], self._ptr_cache["msk"], addr = self.toNumpy()
            nbcells = addr[0]
            dims = addr[5:]
            dims.shape = nbcells, 4
            self._ptr_cache["nbpt"] = dims[:, 0:2]  # number of points, subpoints
            nbval = self._ptr_cache["nbpt"].prod(axis=1)  # nbval by cells (for one component)
            start = np.append([0], nbval.cumsum()[:-1])
            end = start + nbval
            indexes = np.ones((nbcells, nbval.max()), dtype=int) * -1
            for iv in range(nbcells):
                indexes[iv, : nbval[iv]] = np.arange(start[iv], end[iv])
            self._ptr_cache["idx"] = indexes
        return self._ptr_cache

    def __getattr__(self, component):
        """Convenient shortcut to `getComponentValues()`."""
        return self.getComponentValues(component)

    def getComponentValues(self, component: str):
        """Extract the values of a component.

        Args:
            component (str): Component name. Raises ValueError if the component
                does not exist in the field.
        """
        if self._cache["readonly"]:
            self._ptr_cache = None
        icmp = self._cache["cmp"].index(component)
        self._cache["readonly"] = False
        return ComponentOnCells(
            self._cache["val"][:, icmp].copy(),
            ComponentOnCells.Description(
                self._cache["idx"], self._cache["nbpt"], self._cache["msk"][:, icmp]
            ),
        )

    def setComponentValues(self, component: str, cfvalue: ComponentOnCells):
        """Assign the values of a component.

        Args:
            component (str): Component name. Raises ValueError if the component
                does not exist in the field.
            cfvalue (ComponentOnCells): Previously extracted component.
        """
        icmp = self._cache["cmp"].index(component)
        # it directly overwrites '.CESV' vector in place
        self._cache["val"][:, icmp] = cfvalue.expand().getValues()

    def getValues(self, copy=True):
        """
        Returns two numpy arrays containing the field values on specific cells.

        The first array contains the field values while the second one is a mask
        which is `True` if the corresponding value exists, `False` otherwise.
        Array shape is ( number_of_cells_with_components, number_of_components ).

        Args:
            copy (bool): If *True* the data are copied, the is the default.

        Returns:
            ndarray (float): Field values.
            ndarray (bool): Mask for the field values.
        """
        values, mask = self._cache["val"], self._cache["msk"]
        if copy or self._cache["readonly"] is False:
            return values.copy(), mask.copy()

        self._cache["readonly"] = True
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

        values, mask = self._cache["val"], self._cache["msk"]

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
