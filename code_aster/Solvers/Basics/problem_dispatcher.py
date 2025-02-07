# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
from inspect import isabstract

from ...Utilities import logger


class ProblemType(IntFlag):
    """Types of physical problems."""

    Unset = 0
    MecaStat = auto()
    MecaDyna = auto()
    Thermal = auto()


class DispatcherMixin:
    """Mixin class that provides a factory depending on the type of physical problem."""

    @classmethod
    def factory(cls, context):
        """Factory that creates the appropriate object.

        Args:
            context (Context): Context of the problem.

        Returns:
            instance: A new object of the relevant type.
        """
        for kls in cls.__subclasses__():
            logger.debug("candidate: %s", kls)
            if kls.problem_type == context.problem_type:
                if isabstract(kls):
                    return kls.factory(context)
                return kls.builder(context)
        raise TypeError(f"no candidate for {cls=}, type: {context.problem_type}")
