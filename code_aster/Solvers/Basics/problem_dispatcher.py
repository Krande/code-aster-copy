# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

from enum import IntFlag, auto


class ProblemType(IntFlag):
    """Types of physical problems."""

    Unset = 0
    MecaStat = auto()
    MecaDyna = auto()
    Thermal = auto()


class ProblemDispatcher:
    """class that implements the problem dispatcher
    based on the type of physical problem."""

    problem_type = ProblemType.Unset

    @classmethod
    def create(cls, param):
        """factory function that returns the appropriate child class."""
        found = None
        pb_type = cls._getProblemType(param)
        for klass in cls.__subclasses__():
            if klass.problem_type == pb_type:
                found = klass.create(param=param)
        assert found, f"not found: {pb_type.name}"
        found._createPrivate(param)
        return found

    @classmethod
    def _getProblemType(cls, param):
        """Identify the physical problem."""
        if "TYPE_CALCUL" in param:
            pb_type = ProblemType.Thermal
        elif "SCHEMA_TEMPS" in param:
            pb_type = ProblemType.MecaDyna
        else:
            pb_type = ProblemType.MecaStat
        return pb_type

    def _createPrivate(self, param):
        raise NotImplementedError
