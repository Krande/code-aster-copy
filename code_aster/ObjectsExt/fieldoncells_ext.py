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

from libaster import FieldOnCellsReal, FieldOnCellsLong, FieldOnCellsChar8, FieldOnCellsComplex
from ..Objects.Serialization import InternalStateBuilder
from ..Utilities import injector, deprecated


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

    def getValuesWithDescription(self, component, groups=[]):
        """Return the values of a component of the field.

        Arguments:
            component (str): Extracted component.
            groups (list[str], optional): The extraction is limited to the given
                groups of cells.

        Returns:
            tuple(values, description): List of values and description.
            The description provides a tuple with identifiers of
            (cells, points, subpoints).
        """
        mesh = self.getMesh()
        if groups:
            cells = set()
            for grp in groups:
                if mesh.hasGroupOfCells(grp):
                    cells.update(mesh.getCells(grp))
                else:
                    raise ValueError("no {} group of cell".format(grp))
            cells = sorted(cells)
        else:
            cells = mesh.getCells()
        return self.toSimpleFieldOnCells().getValuesWithDescription(cells, component)

    def changePhysicalQuantity(self, physQuantity, map_cmps, fed=None):
        """Return a new field with a new physical quantity and renamed components.

        Arguments:
            physQuantity [str]: name of the new physical quantity
            map_cmps[dict[str, str]]: dict to rename components
            (only renamed component will be keeped)
            fed [FiniteElementDescriptor] : FED used to convert the field if the one present
             in the field is not appropriate (default: None)

        Returns:
            FieldOnCellsReal: field with name physical quantity.
        """

        fcs = self.toSimpleFieldOnCells().changePhysicalQuantity(physQuantity, map_cmps)

        ligrel = self.getDescription()
        if fed:
            ligrel = fed

        return fcs.toFieldOnCells(ligrel)

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
