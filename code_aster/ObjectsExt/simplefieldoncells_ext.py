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

"""
:py:class:`SimpleFieldOnCellsReal`
Simple Fields defined on cells of elements
********************************************************************
"""

try:
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

import warnings

import numpy as np
from libaster import SimpleFieldOnCellsReal

from ..Objects.Serialization import InternalStateBuilder
from ..Utilities import ParaMEDMEM as PMM
from ..Utilities import force_list, injector, logger, no_new_attributes
from ..Utilities import medcoupling as medc
from ..Utilities.mpi_utils import mpi_min, mpi_max, mpi_sum, useHPCMode


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

    Use *in place* operators (`+=`, `*=`...) when possible for performance reasons.

    Args:
        values (ndarray[float]): Values of the component,
            dimension \Sum_i(nbpts_i * nbcmp_i).
        descr (*Description*): Description of the component.
    """

    # WARNING: _idx must be used to access values only on a not restricted component.
    # FIXME '_mask' toujours True
    # '_discr' is not used but interesting in `Description.__repr__`.

    class Description:
        """Internal description of a component.

        Args:
            indexes (ndarray[int]): Indexes in *values* by cells and points on cells.
            discr (ndarray[int]): Local dicretization: number of points and subpoints
                by cell.
            mask (ndarray[bool]): Indicator that tells a value exist, same dimension.
            inner (ndarray[bool], optional): Indicator for inner cells
                (only pass for restricted object).
        """

        _parent = _idx = _mask = _discr = _cells = _nbcells = _nbval = _sign = None
        _inner = None
        __setattr__ = no_new_attributes(object.__setattr__)

        def __init__(self, parent, indexes, discr, mask, inner=None):
            self._parent = parent  # parent SimpleFieldOnCells
            self._idx = indexes
            self._discr = discr
            self._mask = mask
            self._cells = None
            # number of cells in this object before restriction
            self._nbcells = indexes.shape[0]
            # number of values stored for the component before any restriction
            self._nbval = None
            # mask for inner cells (applied on '_cells')
            if inner is not None:
                self._inner = inner
            elif useHPCMode():
                self._inner = np.zeros(self._nbcells, dtype=bool)
                self._inner[parent.getMesh().getInnerCells()] = True
            else:
                self._inner = np.ones(self._nbcells, dtype=bool)
            self._sign = None

        def __sizeof__(self):
            """Return the size of the data in bytes."""
            # for sys.getsizeof(obj)
            return self.sizeof()

        def sizeof(self):
            """Return the size of the data in bytes."""
            return (
                self._idx.nbytes
                + self._discr.nbytes
                + self._mask.nbytes
                + self._inner.nbytes
                + len(list(self._cells if self._cells is not None else [])) * 8
            )

        def setValuesShape(self, shape):
            """Register shape of values if needed."""
            if self._nbval is None:
                self._nbval = shape[0]

        @property
        def sign(self):
            """int: Attribute that holds the signature of the description."""
            if self._sign is None:
                self._sign = (self._idx.ravel() * np.arange(self._idx.size)).sum()
            return self._sign

        def __repr__(self):
            idx = self._idx[self._idx >= 0]
            if len(idx) == 0:
                idx = np.zeros(1, dtype=int)
            lines = [
                f"- indexes shape {self._idx.shape!r}, range {idx.min()} to {idx.max()}, "
                f"{len(idx)} indexes defined, {(self._mask == False).sum()} masked"
            ]
            nbinn = self._inner.sum()
            if self._cells is None:
                lines.append(f"- on all the {self._nbcells} cells, {nbinn} inner")
            else:
                limit = 20
                ndef = len(self._cells)
                lines.append(
                    f"- on {ndef}/{self._nbcells} cells ({nbinn} inner): "
                    + f"{tuple(self._cells[:limit])}"
                    + ("..." if ndef > limit else "")
                )
            uniq, numb = np.unique(self._discr, return_counts=True, axis=0)
            nbv = (uniq.prod(axis=1) * numb).sum()
            discr = [f"{n}x{tuple(d)}" for n, d in zip(numb, uniq)]
            lines.append(
                f"- {len(self._discr)} points/subpoints defs: "
                + ", ".join(discr)
                + f", {nbv} values"
            )
            max = self._idx.max(axis=1)
            undef = max[max < 0]
            nundef = undef.shape[0]
            nb = self._idx.shape[0]
            lines.append(f"- {nb-nundef}/{nb} cells with values")
            return "\n".join(lines)

        def copy(self):
            """Copy of the current description.

            Returns:
                *Decription*: copy of the current description.
            """
            new = __class__(
                self._parent,
                self._idx.copy(),
                self._discr.copy(),
                self._mask.copy(),
                self._inner.copy(),
            )
            new._cells = None if self._cells is None else self._cells.copy()
            new._nbcells = self._nbcells
            new._nbval = self._nbval
            new._sign = self._sign
            return new

        # MPI may be empty
        def restrict(self, cells_ids):
            """Restrict the description on given cells.

            Args:
                cells_ids (list[int]): sublist of the existing cells.

            Returns:
                ndarray[int]: Indexes restricted on the cells list.
            """
            assert self._cells is None, "must be expanded first"
            if isinstance(cells_ids, int):
                cells_ids = force_list(cells_ids)
            self._idx = self._idx[cells_ids]
            idx = self._idx.ravel()
            idx = idx[idx >= 0]
            self._discr = self._discr[cells_ids]
            self._mask = self._mask[idx]
            self._cells = np.array(cells_ids, dtype=int)
            self._inner = self._inner[cells_ids]
            self._sign = None
            return idx

        def masked_values(self, values, undefined=0.0):
            """Returns the vector of masked values.

            Args:
                values (ndarray): Vector of values, restricted or not.
                undefined (float: Value used for undefined components.
            """
            return np.where(self._mask == True, values, undefined)

        def expanded_values(self, values):
            """Expands the vector of values.

            Args:
                values (ndarray): Vector of values, consistent with mask and idx.

            Returns:
                ndarray: Vector of values on all cells (extended with 0.)
            """
            if self._cells is not None:
                values_in = values
                values = np.zeros(self._nbval, dtype=float)
                idx = self._idx[self._idx >= 0]
                values[idx] = values_in
            return values

        def inner_values(self, values):
            """Returns the vector of inner values.

            Args:
                values (ndarray): Vector of values, consistent with mask and idx.

            Returns:
                ndarray: Vector of values on inner cells only (extended with 0.
                    if necessary).
            """
            values = self.expanded_values(self.masked_values(values, undefined=0.0))
            idx = self._idx[self._inner].ravel()
            idx = idx[idx >= 0]
            return values[idx]

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
            return new_values, ComponentOnCells.Description(
                self._parent, new_idx, new_discr, new_mask
            )

    _values = _descr = _sign = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, values, description: Description):
        self.init_from(values, description)

    def init_from(self, values, description: Description):
        """Initialize content."""
        self._values = values
        self._descr = description
        self._descr.setValuesShape(values.shape)
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

    def __sizeof__(self):
        """Return the size of the data in bytes."""
        # for sys.getsizeof(obj)
        return self.sizeof()

    @mpi_sum
    def sizeof(self):
        """Return the size of the data in bytes."""
        return self._values.nbytes + self._descr.sizeof()

    @property
    def innerValues(self):
        """Returns the vector of local values (its shape is inconsistent
        with the internal description)."""
        return self._descr.inner_values(self.values)

    @property
    def values(self):
        """Returns the vector of values (undefined values set to zero)."""
        return self._descr.masked_values(self._values, undefined=0.0)

    @values.setter
    def values(self, values):
        """Set the vector values (components values that are not expected
        are forced to 0.).

        Args:
            values (ndarray[float]): New values.
        """
        assert values.shape == self._values.shape
        self._values = self._descr.masked_values(values)

    def getValuesByCells(self):
        """Returns the values by cells, expanded by 0.0 where undefined."""
        if self._cells is not None:
            return self.expand().getValuesByCells()[self._cells]
        values = np.append(self.values, [0.0])
        return np.take(values, self._idx)

    def getCells(self):
        """Return the cells id."""
        if self._cells is None:
            return np.arange(self._descr._nbcells)
        return self._cells

    def getNumberOfCells(self):
        """Return the total number of cells (even if the current component is
        currently restricted).

        Returns:
            int: Number of cells the field is defined on.
        """
        return self._descr._nbcells

    @property
    def size(self):
        """Return the number of inner values of the vector."""
        return len(self.innerValues)

    def __repr__(self):
        lines = repr(self._values).splitlines()
        if len(lines) > 6:
            lines = lines[:3] + ["..."] + lines[-3:]
        lines.insert(0, f"- {len(self._values)} values in vector:")
        lines.append(repr(self._descr))
        return "\n".join(lines)

    def plot_descr(self, filename=None):
        """Plot the description of values by cells.

        A matplotlib GUI backend should be enabled, for example with:
        ``matplotlib.use("TkAgg")``.

        Args:
            filename (str|Path, optional): Save the image to this file if provided.
                Otherwise plot interactively the figure.
        """
        if not HAS_MATPLOTLIB:
            warnings.warn("matplotlib is not available")
            return
        idx = np.where(self._idx >= 0, 1, 0)
        undef = np.where(self._descr._mask[self._idx] == False, -2, 0)
        idx = np.where(undef == 0, idx, undef * idx)
        fig, ax = plt.subplots()
        ax.grid(True, linestyle="--")
        ax.set_xlabel("numbers of points")
        ax.set_ylabel("cells id")
        cmap = ListedColormap([(0.6, 0.6, 0.6), "white", (0.0, 0.4, 0.8)], "indexed")
        plt.imshow(idx, interpolation="nearest", cmap=cmap, aspect="auto", vmin=-1, vmax=1)
        if filename:
            plt.savefig(filename)
        else:
            plt.show()

    @property
    def sign(self):
        """int: Attribute that holds the signature of the component."""
        if self._sign is None:
            self._sign = self._descr.sign + len(self._values)
        return self._sign

    def _check_consistency(self, other):
        """Check consistency between the both signatures. A float is returned as is."""
        if isinstance(other, ComponentOnCells):
            if self.sign != other.sign:
                logger.warning("lhs: %s", repr(self))
                logger.warning("rhs: %s", repr(other))
                raise IndexError("inconsistent shapes")
            other = other._values
        return other

    def apply(self, func):
        """Apply a function of the values."""
        if not isinstance(func, np.ufunc):
            func = np.vectorize(func)
        return ComponentOnCells(func(self._values), self._descr.copy())

    def applyOnArray(self, func):
        """Apply a function of the array of values."""
        return ComponentOnCells(func(self._values), self._descr.copy())

    def __iadd__(self, other):
        other = self._check_consistency(other)
        self._values += other
        return self

    def __add__(self, other):
        other = self._check_consistency(other)
        return ComponentOnCells(self._values + other, self._descr.copy())

    def __radd__(self, other):
        return self + other

    def __isub__(self, other):
        other = self._check_consistency(other)
        self._values -= other
        return self

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __neg__(self):
        return ComponentOnCells(-self._values, self._descr.copy())

    def __imul__(self, other):
        other = self._check_consistency(other)
        self._values *= other
        return self

    def __mul__(self, other):
        other = self._check_consistency(other)
        return ComponentOnCells(self._values * other, self._descr.copy())

    def __rmul__(self, other):
        return self * other

    def __itruediv__(self, other):
        other = self._check_consistency(other)
        other = self._values(other, 1.0)
        if np.any(other == 0.0):
            raise ZeroDivisionError()
        self._values /= other
        return self

    def __truediv__(self, other):
        other = self._check_consistency(other)
        other = self._descr.masked_values(other, 1.0)
        if np.any(other == 0.0):
            raise ZeroDivisionError()
        return ComponentOnCells(self._values / other, self._descr.copy())

    def __rtruediv__(self, other):
        values = self._descr.masked_values(self._values, 1.0)
        if np.any(values == 0.0):
            raise ZeroDivisionError()
        return ComponentOnCells(1.0 / values * other, self._descr.copy())

    def __abs__(self):
        return ComponentOnCells(np.abs(self._values), self._descr.copy())

    def __ipow__(self, other):
        other = self._check_consistency(other)
        self._values **= other
        return self

    def __pow__(self, other):
        other = self._check_consistency(other)
        return ComponentOnCells(np.power(self._values, other), self._descr.copy())

    @mpi_min
    def min(self):
        """Return the minimum value.

        Returns:
            float: minimum of the values.
        """
        return self.innerValues.min()

    @mpi_max
    def max(self):
        """Return the maximum value.

        Returns:
            float: maximum of the values.
        """
        return self.innerValues.max()

    @mpi_sum
    def sum(self):
        """Return the sum of the values.

        Returns:
            float: sum of the values.
        """
        return self.innerValues.sum()

    def mean(self):
        """Return the mean value.

        Returns:
            float: mean of the values.
        """
        return self.innerValues.mean()

    def minimum(self, other):
        """Compare the component with another one and keep the element-wise minima.

        Args:
            other (*ComponentOnCells*): Another component.
        """
        other = self._check_consistency(other)
        self._values = np.minimum(self._values, other)

    def maximum(self, other):
        """Compare the component with another one and keep the element-wise maxima.

        Args:
            other (*ComponentOnCells*): Another component.
        """
        other = self._check_consistency(other)
        self._values = np.maximum(self._values, other)

    # @DebugArgs.pickle_on_error
    def restrict(self, cells_ids):
        """Restrict the component on given cells.

        Args:
            cells_ids (list[int]): Cells ids to be used.
        """
        if self._cells is not None:
            exp = self.expand()
            self.init_from(exp._values, exp._descr)
        idx = self._descr.restrict(cells_ids)
        self._values = self._values[idx]
        self._sign = None

    def filterByValues(self, mini, maxi, strict_mini=True, strict_maxi=True):
        """Returns the cells ids where at least one value of the cell is in
        the interval.

        Args:
            mini (float): minimum value.
            maxi (float): maximum value.
            strict_mini (bool): Use strict comparison (the default) for the minimum.
            strict_maxi (bool): Use strict comparison (the default) for the maximum.

        Returns:
            list[int]: Cells ids.
        """
        if self._cells is not None:
            cells = self.expand().filterByValues(mini, maxi, strict_mini, strict_maxi)
            return list(set(self._cells).intersection(cells))
        if strict_mini:
            filter = mini < self._values
        else:
            filter = mini <= self._values
        if strict_maxi:
            filter &= self._values < maxi
        else:
            filter &= self._values <= maxi

        filter = np.append(filter, [False])  # not defined (-1): not in interval...
        idx = np.where(self._descr._mask[self._idx] == False, -1, self._idx)  # or masked
        idx_filt = filter[idx]
        cells_filt = idx_filt.max(axis=1)
        cells_filt &= self._descr._inner
        cells = list(np.arange(self.getNumberOfCells(), dtype=int)[cells_filt])
        return cells

    # @DebugArgs.pickle_on_error
    def expand(self):
        """Expand the component by adding 0. where it not defined.

        Returns:
            *ComponentOnCells*: expanded component.
        """
        if self._cells is None:
            return self.copy()
        return ComponentOnCells(*self._descr.expand(self._values))

    # @DebugArgs.pickle_on_error
    def onSupportOf(self, other, strict=False):
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
        if self._cells is None and other._cells is not None:
            changed = self.onSupportOf(other.expand())
            changed.restrict(other._cells)
            return changed
        if self._cells is not None and other._cells is None:
            return self.expand().onSupportOf(other)
        if self._cells is not None:
            if np.any(self._cells != other._cells):
                raise IndexError("components must be defined on the same cells")
            return self.expand().onSupportOf(other)
        assert self._descr._nbcells == other._descr._nbcells, (
            self._descr._nbcells == other._descr._nbcells
        )
        assert self._cells is None and other._cells is None, (self._cells, other._cells)
        assert self._idx.shape == other._idx.shape
        new = other.copy()
        bools = other._idx >= 0
        idx = other._idx[bools]
        keep = self._idx[bools]
        if strict and np.any(keep < 0):
            raise IndexError("no value on all the support")
        prol = np.append(self._values, [0.0])
        new._values[idx] = np.take(prol, keep)
        new._descr._idx = other._idx.copy()
        self._sign = None
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
                self, self._cache["idx"], self._cache["nbpt"], self._cache["msk"][:, icmp]
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
        self._cache["val"][:, icmp] = cfvalue.expand().values

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
