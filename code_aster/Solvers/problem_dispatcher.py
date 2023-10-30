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


class ProblemType:
    """Types of physical problems."""

    Unset = 0x00
    Static = 0x01
    MecaDyna = 0x02
    Thermal = 0x04


class ProblemDispatcher:
    """xxxx xxxx xxxx xxxx."""

    problem_type = ProblemType.Unset

    @classmethod
    def create(cls, param):
        """xxxx xxxx xxxx xxxx."""
        found = None
        ProblemType = cls._getProblemType(param)
        for klass in cls.__subclasses__():
            if klass.problem_type == ProblemType:
                found = klass.create(param=param)
        if found is None:
            raise RuntimeError
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
            pb_type = ProblemType.Static
        return pb_type

    def _createPrivate(self, param):
        """xxxx xxxx xxxx xxxx."""
        raise NotImplementedError
