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

"""
:py:class:`FieldOnCellsReal` --- Fields defined on nodes of elements
********************************************************************

The *Field On Nodes* object exists for real numbers (:py:class:`FieldOnCellsReal`),
integers (:py:class:`FieldOnCellsLong`),
strings (:py:class:`FieldOnCellsChar8`) and
complex numbers (:py:class:`FieldOnCellsComplex`).
"""

import os.path as osp
import subprocess

from libaster import FieldOnCellsReal, FieldOnCellsLong, FieldOnCellsChar8, FieldOnCellsComplex
from ..Objects.Serialization import InternalStateBuilder
from ..Utilities import injector, force_list
from ..Utilities import MPI, ExecutionParameter, SharedTmpdir


def _getValuesWithDescription_numpy(sfield, cmps_arg, groups):
    """Pure-Python fallback for getValuesWithDescription using toNumpy().

    When the C++ loop returns empty despite valid data (Windows bug),
    this function extracts the same information from the raw Jeveux
    arrays exposed by toNumpy().

    Returns the same format as getValuesWithDescription:
        (values, (cells, cmps, points, subpoints))
    """
    ncmps = sfield.getNumberOfComponents()
    if ncmps == 0:
        return ([], ([], [], [], []))

    all_cmp_names = sfield.getComponents()
    cmp_name2idx = {name: i for i, name in enumerate(all_cmp_names)}

    # Determine which components to extract
    if cmps_arg:
        cmp_indices = []
        cmp_names_out = []
        for c in cmps_arg:
            if c in cmp_name2idx:
                cmp_indices.append(cmp_name2idx[c])
                cmp_names_out.append(c)
    else:
        cmp_indices = list(range(ncmps))
        cmp_names_out = list(all_cmp_names)

    if not cmp_indices:
        return ([], ([], [], [], []))

    # Determine which cells to iterate
    mesh = sfield.getMesh()
    grp_list = force_list(groups) if groups else []
    if grp_list and all(isinstance(g, (int,)) for g in grp_list):
        cells = grp_list
    else:
        cells = list(mesh.getCells(grp_list) if grp_list else mesh.getCells())

    # Get raw numpy arrays: values(n_rows, ncmps), mask(n_rows, ncmps), size_arr(1D)
    vals_arr, mask_arr, size_arr = sfield.toNumpy()

    values = []
    v_cells = []
    v_cmps = []
    points = []
    subpoints = []

    for cell in cells:
        npt = int(size_arr[4 + 4 * cell + 1])
        nspt = int(size_arr[4 + 4 * cell + 2])
        cell_ncmp = int(size_arr[4 + 4 * cell + 3])
        decal = int(size_arr[4 + 4 * cell + 4])

        if npt == 0 or nspt == 0 or cell_ncmp == 0:
            continue

        # decal is always a multiple of ncmps (cell_ncmp == ncmps for allocated cells)
        base_row = decal // ncmps

        for icmp, cmp_name in zip(cmp_indices, cmp_names_out):
            if icmp >= cell_ncmp:
                continue
            for ipt in range(npt):
                for ispt in range(nspt):
                    row = base_row + ipt * nspt + ispt
                    col = icmp
                    if row < mask_arr.shape[0] and mask_arr[row, col]:
                        values.append(float(vals_arr[row, col]))
                        v_cells.append(cell)
                        v_cmps.append(cmp_name)
                        points.append(ipt)
                        subpoints.append(ispt)

    return (values, (v_cells, v_cmps, points, subpoints))


class FieldOnCellsStateBuilder(InternalStateBuilder):
    """Class that returns the internal state of a *FieldOnCells*."""

    def save(self, field):
        """Return the internal state of a *Result* to be pickled.

        Arguments:
            field (*FieldOnCells*): The *FieldOnCells* object to be pickled.

        Returns:
            *InternalStateBuilder*: The internal state itself.
        """
        super().save(field)
        self._st["fed"] = field.getDescription()
        return self

    def restore(self, field):
        """Restore the *DataStructure* content from the previously saved internal
        state.

        Arguments:
            field (*DataStructure*): The *DataStructure* object to be restored.
        """
        super().restore(field)
        if self._st["fed"]:
            field.setDescription(self._st["fed"])


@injector(FieldOnCellsReal)
class ExtendedFieldOnCellsReal:
    cata_sdj = "SD.sd_champ.sd_cham_elem_class"
    internalStateBuilder = FieldOnCellsStateBuilder

    def getValuesWithDescription(self, components=[], groups=[]):
        """Return the values of a component of the field.

        Arguments:
            components (list[str]): Extracted components.
            groups (list[str], optional): The extraction is limited to the given
                groups of cells.

        Returns:
            tuple(values, description): List of values and description.
            The description provides a tuple with identifiers of
            (cells, points, subpoints).
        """

        import sys
        sfield = self.toSimpleFieldOnCells()
        cmps_arg = force_list(components)
        grp_list = force_list(groups) if groups else []
        # When groups contains integers, they are cell indices (not group
        # names).  The C++ getValuesWithDescription expects group name
        # strings, so we must use the numpy fallback for integer indices.
        if grp_list and all(isinstance(g, (int,)) for g in grp_list):
            return _getValuesWithDescription_numpy(sfield, cmps_arg, groups)
        result = sfield.getValuesWithDescription(cmps_arg, groups)
        # Workaround: on Windows, the C++ getValuesWithDescription loop can
        # return empty despite valid data.  Fall back to a pure-Python
        # extraction using the raw numpy arrays from toNumpy().
        if len(result[0]) == 0 and sys.platform == "win32":
            result = _getValuesWithDescription_numpy(sfield, cmps_arg, groups)
        return result

    def plot(self, command="gmsh", local=False, split=False):
        """Plot the field.

        Arguments:
            command (str): Program to be executed to plot the field.
            local (bool): Print in separate files if *True*. Otherwise an unique file is used.
            split (bool): Display the field on each subdomain separately if *True*. Otherwise the global field is displayed.
        """
        comm = MPI.ASTER_COMM_WORLD
        mesh = self.getMesh()
        opt = "Mesh.VolumeEdges = 0;Mesh.VolumeFaces=0;Mesh.SurfaceEdges=0;Mesh.SurfaceFaces=0;View[0].ShowElement = 1;"
        if (local and mesh.isParallel()) or (split and mesh.isParallel()):
            with SharedTmpdir("plot") as tmpdir:
                filename = osp.join(tmpdir.path, f"field_{comm.rank}.med")
                self.printMedFile(filename, local=True)
                comm.Barrier()
                if comm.rank == 0:
                    if split:
                        for i in range(comm.size):
                            ff = osp.join(tmpdir.path, f"field_{i}.med")
                            subprocess.run(
                                [
                                    ExecutionParameter().get_option(f"prog:{command}"),
                                    "-string",
                                    opt,
                                    ff,
                                ]
                            )
                    else:
                        files = [osp.join(tmpdir.path, f"field_{i}.med") for i in range(comm.size)]
                        subprocess.run(
                            [ExecutionParameter().get_option(f"prog:{command}"), "-string", opt]
                            + files
                        )
        else:
            with SharedTmpdir("plot") as tmpdir:
                filename = osp.join(tmpdir.path, "field.med")
                self.printMedFile(filename, local=False)
                if comm.rank == 0:
                    subprocess.run(
                        [
                            ExecutionParameter().get_option(f"prog:{command}"),
                            "-string",
                            opt,
                            filename,
                        ]
                    )
        print("waiting for all plotting processes...")
        comm.Barrier()

    def asPhysicalQuantity(self, physQuantity, map_cmps, fed=None):
        """Return a new field with a new physical quantity and renamed components.

        Arguments:
            physQuantity [str]: name of the new physical quantity
            map_cmps [dict[str, str]]: dict to rename components
                (only renamed components will be keeped)
            fed [FiniteElementDescriptor] : FED used to convert the field if
                the one present in the field is not appropriate (default: None)

        Returns:
            FieldOnCellsReal: field with name physical quantity.
        """
        fcs = self.toSimpleFieldOnCells().asPhysicalQuantity(physQuantity, map_cmps)
        ligrel = self.getDescription()
        if fed:
            ligrel = fed
        return fcs.toFieldOnCells(ligrel)


@injector(FieldOnCellsLong)
class ExtendedFieldOnCellsLong:
    cata_sdj = "SD.sd_champ.sd_cham_elem_class"
    internalStateBuilder = FieldOnCellsStateBuilder


@injector(FieldOnCellsChar8)
class ExtendedFieldOnCellsChar8:
    cata_sdj = "SD.sd_champ.sd_cham_elem_class"
    internalStateBuilder = FieldOnCellsStateBuilder


@injector(FieldOnCellsComplex)
class ExtendedFieldOnCellsComplex:
    cata_sdj = "SD.sd_champ.sd_cham_elem_class"
    internalStateBuilder = FieldOnCellsStateBuilder
