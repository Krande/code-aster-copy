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

from ..Supervis import ExecuteCommand


class ModiRepere(ExecuteCommand):
    """Command MODI_REPERE
    """
    command_name = "MODI_REPERE"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        if "reuse" in keywords:
            self._result = keywords["reuse"]
        elif "RESULTAT" in keywords:
            if keywords["RESULTAT"].getType() == "EVOL_ELAS":
                self._result = type(keywords["RESULTAT"])()
            else:
                self._result = type(keywords["RESULTAT"])()
        else:
            self._result = type(keywords["CHAM_GD"])()

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        modele = None
        try:
            modele = keywords["RESULTAT"].getModel()
            if modele is not None:
                self._result.setModel(modele)
                self._result.build()
        except:
            pass

    def add_dependencies(self, keywords):
        """Register input *DataStructure* objects as dependencies.

        Arguments:
            keywords (dict): User's keywords.
        """
        super().add_dependencies(keywords)
        self.remove_dependencies(keywords, "RESULTAT")
        self.remove_dependencies(keywords, "CHAM_GD")


MODI_REPERE = ModiRepere.run
