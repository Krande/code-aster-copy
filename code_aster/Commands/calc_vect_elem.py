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

from ..Objects import (ElementaryVectorDisplacementReal,
                       ElementaryVectorPressureComplex,
                       ElementaryVectorTemperatureReal)
from ..Supervis import ExecuteCommand
from ..Utilities import force_list


class ComputeElementaryVector(ExecuteCommand):

    """Command that creates elementary vectors."""
    command_name = "CALC_VECT_ELEM"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        if keywords['OPTION'] == "CHAR_MECA": self._result = ElementaryVectorDisplacementReal()
        elif keywords['OPTION'] == "CHAR_THER": self._result = ElementaryVectorTemperatureReal()
        elif keywords['OPTION'] == "CHAR_ACOU": self._result = ElementaryVectorPressureComplex()
        else: raise NotImplementedError("Must be implemented")

    def post_exec(self, keywords):
        """Store references to ElementaryVector objects.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords, changed
                in place.
        """
        if "MODELE" in keywords:
            self._result.setModel(keywords["MODELE"])

        if "CARA_ELEM" in keywords:
            self._result.setElementaryCharacteristics(keywords["CARA_ELEM"])

        if "CHAM_MATER" in keywords:
            self._result.setMaterialField(keywords["CHAM_MATER"])

        if "CHARGE" in keywords:
            loads = force_list(keywords["CHARGE"])
            self._result.setModel(loads[0].getModel())

        self._result.build()

CALC_VECT_ELEM = ComputeElementaryVector.run
