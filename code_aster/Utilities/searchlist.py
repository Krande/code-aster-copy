# coding: utf-8
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

# person_in_charge: francesco.bettonte@edf.fr

"""
:py:mod:`searchlist` --- General purpose utilities for search lines
***************************************************************

This modules gives some basic utilities for search lines.
"""

import numpy as np


class SearchList:
    @property
    def values(self):
        """Attribute that holds the list values.

        Returns:
            list[float|int]: List of values
        """
        return self._values

    @property
    def precision(self):
        """Attribute that holds the search precision.

        Returns:
            float: The search precision
        """
        return self._precision

    @precision.setter
    def precision(self, value):
        assert value > np.finfo(float).eps
        self._precision = value

    @property
    def criterion(self):
        """Attribute that holds the search criterion.

        Returns:
            str: The search criterion, RELATIF | ABSOLU
        """
        return self._criterion

    @criterion.setter
    def criterion(self, value):
        expected = ("RELATIF", "ABSOLU")
        if value not in expected:
            msg = "Criterion '{0}' is not in {1}".format(value, expected)
            raise ValueError(msg)
        self._criterion = value

    def __repr__(self):
        return "[{0}]".format(", ".join(map(str, self.values)))

    def __iter__(self):
        yield from self.values

    def __init__(self, values, precision=1.0e-6, criterion="RELATIF"):
        """Initialization of the Search List

        Arguments:
            values (list[float|int]): The list of values for search
            precision (float): The search precision
            criterion (str): The search criterion ( ABSOLU|RELATIF )
        """

        assert all(isinstance(i, (int, float)) for i in values)
        self._values = np.array(values)
        self.precision = precision
        self.criterion = criterion

    def _search_candidates(self, value):
        """
        Search indexes of value matching the defined criterion
        """
        assert isinstance(value, (int, float))

        if self.criterion == "RELATIF":
            min_v = value * (1 - self.precision)
            max_v = value * (1 + self.precision)
        else:
            min_v = value - self.precision
            max_v = value + self.precision

        min_v, max_v = sorted((min_v, max_v))
        return np.flatnonzero(np.logical_and(self.values >= min_v, self.values <= max_v))

    def __contains__(self, value):
        idx = self._search_candidates(value)
        return len(idx) > 0

    def unique(self, value):
        """
        Return True if value is unique
        """
        idx = self._search_candidates(value)
        return len(idx) == 1

    def index(self, value):
        """
        Return index of value if value is unique
        """
        idx = self._search_candidates(value)
        if len(idx) == 0:
            msg = "{0} is not in list".format(value)
            raise ValueError(msg)
        elif len(idx) > 1:
            msg = "{0} is not unique in list".format(value)
            raise IndexError(msg)
        else:
            return idx[0]
