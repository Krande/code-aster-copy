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

from ..Objects import ModeResult
from ..Supervis import ExecuteCommand


class ModalSeismicCombination(ExecuteCommand):
    """Command that creates the :class:`~code_aster.Objects.NotImplemented`"""

    command_name = "COMB_SISM_MODAL"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        self._result = ModeResult()

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        feds = []
        for option in keywords["OPTION"]:
            if option.endswith("_ELGA") or option.endswith("_ELNO"):
                feds = keywords["MODE_MECA"].getFiniteElementDescriptors()
                break

        dofNum = keywords["MODE_MECA"].getDOFNumbering()
        if dofNum is not None:
            self._result.setDOFNumbering(dofNum)
        model = keywords["MODE_MECA"].getModel()
        if model:
            self._result.setModel(model)
        else:
            self._result.setMesh(keywords["MODE_MECA"].getMesh())
        self._result.build(feds=feds)


COMB_SISM_MODAL = ModalSeismicCombination.run
