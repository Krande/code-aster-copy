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

# person_in_charge: mathieu.courtois@edf.fr
"""
:py:class:`FieldOnCellsReal` --- Fields defined on nodes of elements
********************************************************************

The *Field On Nodes* object exists for real numbers (:py:class:`FieldOnCellsReal`),
integers (:py:class:`FieldOnCellsLong`),
strings (:py:class:`FieldOnCellsChar8`) and
complex numbers (:py:class:`FieldOnCellsComplex`).
"""

import numpy
import os.path as osp
import subprocess

from libaster import FieldOnCellsReal, FieldOnCellsLong, FieldOnCellsChar8, FieldOnCellsComplex
from ..Objects.Serialization import InternalStateBuilder
from ..Utilities import injector, deprecated, force_list
from ..Utilities import MPI, ExecutionParameter, SharedTmpdir


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
        result = sfield.getValuesWithDescription(cmps_arg, groups)
        # Temporary diagnostic: log when getValuesWithDescription returns empty
        # to help debug Windows-specific empty results (see win-support branch).
        if len(result[0]) == 0 and sys.platform == "win32":
            ncells = sfield.getNumberOfCells()
            ncmps = sfield.getNumberOfComponents()
            cmp_names = sfield.getComponents() if ncmps > 0 else []
            maxpt = sfield.getMaxNumberOfPoints()
            maxspt = sfield.getMaxNumberOfSubPoints()
            # Check mesh cells resolution
            mesh = sfield.getMesh()
            grp_list = force_list(groups) if groups else []
            resolved_cells = mesh.getCells(grp_list) if grp_list else mesh.getCells()
            # Sample first few cells for point/subpoint/hasValue info
            sample_info = []
            for c in resolved_cells[:3]:
                try:
                    npt = sfield.getNumberOfPointsOfCell(c)
                    nspt = sfield.getNumberOfSubPointsOfCell(c)
                    hv = sfield.hasValue(c, 0, 0, 0) if npt > 0 and nspt > 0 and ncmps > 0 else None
                    sample_info.append(f"cell{c}:npt={npt},nspt={nspt},hv={hv}")
                except Exception as e:
                    sample_info.append(f"cell{c}:err={e}")
            print(
                f"[DIAG] getValuesWithDescription returned empty: "
                f"ncells={ncells}, ncmps={ncmps}, maxpt={maxpt}, maxspt={maxspt}, "
                f"cmps_arg={cmps_arg}, groups={groups}, "
                f"field_cmps={cmp_names[:8]}, "
                f"resolved_cells_count={len(resolved_cells)}, "
                f"samples=[{'; '.join(sample_info)}]",
                flush=True,
            )
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

    # <16.4.0
    @deprecated(case=1, help="Use 'asLocalization()' instead, it returns a new field.")
    def changeLocalization(self, loc):
        return self.asLocalization(loc)

    # <16.4.0
    @deprecated(case=4, help="Use 'getValuesWithDescription()' instead")
    def EXTR_COMP(self, comp, lgma=[], topo=0):
        """Deprecated: Use 'getValuesWithDescription()' instead.

        Examples:

        .. code-block:: python

            # previously:
            extrcmp = chamele.EXTR_COMP(cmp, groups, 1)
            values = extrcmp.valeurs
            cells = extrcmp.maille
            points = extrcmp.point
            subpoints = extrcmp.sous_point
            # replaced by:
            values, (cells, points, subpoints) = chamele.getValuesWithDescription(cmp, groups)

            # previously:
            extrcmp = chamele.EXTR_COMP(cmp, groups, 0)
            values = extrcmp.valeurs
            # replaced by:
            values, _  = chamele.getValuesWithDescription(cmp, groups)
        """


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
