# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

import aster
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
            field (*DataStructure*): The *DataStructure* object to be pickled.
        """
        super().restore(field)
        if self._st["fed"]:
            field.setDescription(self._st["fed"])


@injector(FieldOnCellsReal)
class ExtendedFieldOnCellsReal:
    cata_sdj = "SD.sd_champ.sd_cham_elem_class"
    internalStateBuilder = FieldOnCellsStateBuilder

    def getValuesWithDescription(self, comp, lgma=[]):
        """Retourne les valeurs de la composante comp du champ sur la liste
        de groupes de mailles lgma avec avec la description.
        Si lgma est une liste vide, c'est equivalent a un TOUT='OUI'
        """
        mesh = self.getMesh()
        if lgma:
            cells = set()
            for grMa in lgma:
                if mesh.hasGroupOfCells(grMa):
                    cells.update(mesh.getCells(grMa))
                else:
                    raise ValueError("no {} group of cell".format(grMa))
            cells = sorted(cells)
        else:
            cells = mesh.getCells()

        return self.toSimpleFieldOnCells().getValuesWithDescription(cells, comp)

    @deprecated
    def EXTR_COMP(self, comp, lgma=[], topo=0):
        raise Exception(
            """EXTR_COMP has been removed, use getValuesWithDescription instead
        Ex1:
            extrcmp = chamele.EXTR_COMP(cmp, groups)
            values = extrcmp.valeurs
            cells = extrcmp.maille
            points = extrcmp.point
            subpoints = extrcmp.sous_point
        =>
            values, (cells, points, subpoints) = chamele.getValuesWithDescription(cmp, groups)
        Ex2:
            extrcmp = chamele.EXTR_COMP(cmp, groups)
            values = extrcmp.valeurs
        =>
            values, _  = chamele.getValuesWithDescription(cmp, groups)
        """
        )


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
