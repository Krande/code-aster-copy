# coding: utf-8

# Copyright (C) 1991 - 2020  EDF R&D                www.code-aster.org
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

from ..Objects import FieldOnNodesReal
from ..Supervis import ExecuteCommand


class SolveLinearSystem(ExecuteCommand):
    """Command that solves :class:`~code_aster.Objects.AssemblyMatrix`."""
    command_name = "RESOUDRE"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        self._result = FieldOnNodesReal()

    def add_dependencies(self, keywords):
        """Register input *DataStructure* objects as dependencies.
        Add dependency to DOFNumbering instead of all Matrix object
        Arguments:
            keywords (dict): User's keywords.
        """
        super().add_dependencies(keywords)

        Matr = keywords.get("MATR")
        Nume = Matr.getDOFNumbering()
        self._result.addDependency(Nume)

        self.remove_dependencies(keywords, "MATR")

RESOUDRE = SolveLinearSystem.run
