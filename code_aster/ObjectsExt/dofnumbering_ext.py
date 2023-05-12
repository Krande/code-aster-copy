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
:py:class:`DOFNumbering` --- DOFNumbering definition
*****************************************
"""

from libaster import DOFNumbering

from ..Utilities import injector
import functools


@injector(DOFNumbering)
class ExtendedDOFNumbering:
    cata_sdj = "SD.sd_nume_ddl.sd_nume_ddl"

    def __getinitargs__(self):
        """Returns the argument required to reinitialize a
        DOFNumbering object during unpickling.
        """
        return (self.getName(), self.getEquationNumbering(), self.getModel())

    @property
    @functools.lru_cache()
    def __Components2Rows(self, local=True):
        """Build the dictionary from the components to the rows."""
        ndofs = self.getNumberOfDOFs()
        dict_dof = {}
        for row in range(ndofs):
            component = self.getComponentFromDOF(int(row))
            dict_dof.setdefault(component, []).append(row)
        return dict_dof

    def getRowsAssociatedToComponent(self, component: str, local=True):
        """Return the rows associated to the input component.

        Arguments:
            component (str): the name of the component (aka degree of freedom)
        """
        available_components = self.getComponents()
        if component not in available_components:
            raise ValueError(f"Component {component} is not in {available_components}")
        return self.__Components2Rows[component]

    def getDictComponentsToRows(self, local=True):
        """Return the dictionary with the available components as keys and the rows as values."""
        return self.__Components2Rows
