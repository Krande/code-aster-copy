# coding: utf-8

# Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

from ..Objects import Model
from ..Supervis import ExecuteCommand


class ModelAssignment(ExecuteCommand):
    """Command that creates the :class:`~code_aster.Objects.Model` by assigning
    finite elements on a :class:`~code_aster.Objects.Mesh`."""
    command_name = "AFFE_MODELE"

    def adapt_syntax(self, keywords):
        """Hook to adapt syntax *after* syntax checking.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords, changed
                in place.
        """
        if keywords["MAILLAGE"].isParallel():
            keywords["DISTRIBUTION"] = {'METHODE': 'CENTRALISE'}

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        self._result = Model(keywords["MAILLAGE"])

        # set model in FED (because of circular reference)
        FED = self._result.getFiniteElementDescriptor()
        FED.setModel(self._result)

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """

        FED = self._result.getFiniteElementDescriptor()
        FED.build()

    def add_dependencies(self, keywords):
        """Register input *DataStructure* objects as dependencies.

        Arguments:
            keywords (dict): User's keywords.
        """

        # Add no dependencies since everything is mesh is already added


AFFE_MODELE = ModelAssignment.run
