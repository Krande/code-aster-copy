# coding: utf-8

# Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

from ..Objects import FieldOnNodesReal
from ..Supervis import ExecuteCommand


class SismicLoading(ExecuteCommand):
    """ """

    command_name = "CALC_CHAR_SEISME"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """

        self._result = FieldOnNodesReal()

    def post_exec(self, keywords):
        """Build.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """

        nume_ddl = keywords["MATR_MASS"].getDOFNumbering()
        self._result.setDescription(nume_ddl.getEquationNumbering())


CALC_CHAR_SEISME = SismicLoading.run
