# coding: utf-8

# Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
#
# This file is part of Code_Aster.
#
# Code_Aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Code_Aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.

# person_in_charge: nicolas.sellenet@edf.fr

from ..Objects import DOFNumbering, ParallelDOFNumbering
from ..Supervis import ExecuteCommand
from ..Utilities import logger


class NumberingCreation(ExecuteCommand):
    """Command that creates a :class:`~code_aster.Objects.DOFNumbering`."""
    command_name = "NUME_DDL"

    def meshIsParallel(self, keywords):
        """ The mesh is a ParallelMesh ?

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """

        model = keywords.get('MODELE')
        if model is not None:
            return model.getMesh().isParallel()
        else:
            matr = keywords.get('MATR_RIGI')[0]
            return matr.getMesh().isParallel()


    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """

        if self.meshIsParallel(keywords):
            self._result = ParallelDOFNumbering()
        else:
            self._result = DOFNumbering()


    def exec_(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """

        model = keywords.get('MODELE')
        if model is not None:
            self._result.setModel(model)
            charge = keywords.get("CHARGE")
            if charge is not None:
                for curLoad in charge:
                    self._result.addLoad(curLoad)
        else:
            matrRigi = keywords['MATR_RIGI']
            for matr in matrRigi:
                self._result.setElementaryMatrix(matr)

        self._result.computeNumbering()


NUME_DDL = NumberingCreation.run
