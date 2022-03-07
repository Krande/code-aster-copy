# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
:py:class:`DiscreteComputation` --- DiscreteComputation object
******************************************************
"""

from libaster import DiscreteComputation

from ..Utilities import injector


@injector(DiscreteComputation)
class ExtendedDiscreteComputation:

    def __getinitargs__(self):
        """Returns the argument required to reinitialize a MaterialField
        object during unpickling.
        """
        return (self.getPhysicalProblem(), )

    def elasticStiffnessMatrix(self, time=0.0, fourierMode=0, groupsOfCells=None):
        """Return the elementary matices for elastic Stiffness matrix

      Arguments:
            time (float): current time (default = 0.0)
            fourierMode (int): Fourier mode (default = 0)
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it is None, the full model is used (default = None)

      Returns:
            ElementaryMatrix: elementary elastic Stiffness matrices
        """
        groups = groupsOfCells
        if groups is None:
            groups = []

        return self.elasticStiffnessMatrix_(time, fourierMode, groups)
