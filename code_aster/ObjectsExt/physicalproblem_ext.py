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
:py:class:`PhysicalProblem` --- PhysicalProblem object
******************************************************
"""

from ..Objects import PhysicalProblem, DirichletBC
from ..Utilities import injector


@injector(PhysicalProblem)
class ExtendedPhysicalProblem:

    def __getinitargs__(self):
        """Returns the argument required to reinitialize a MaterialField
        object during unpickling.
        """
        return (self.getModel(), self.getMaterialField(),
            self.getElementaryCharacteristics())

    def addLoadFromDict(self, dictLoad):
        """Add a load from a dict - to remove quickly

        Arguments:
            a dict with key "CHARGE" for the load and
            optionally "FONC_MULT" for a load function
        """
        charge = dictLoad["CHARGE"]
        if "FONC_MULT" in dictLoad:
            if isinstance(charge, DirichletBC):
                self.addDirichletBC(charge, dictLoad["FONC_MULT"])
            else:
                self.addLoad(charge, dictLoad["FONC_MULT"])
        else:
            if isinstance(charge, DirichletBC):
                self.addDirichletBC(charge)
            else:
                self.addLoad(charge)

        if not hasattr(self, 'allLoadsDict'):
            self.allLoadsDict = []
        self.allLoadsDict.append(dictLoad)
