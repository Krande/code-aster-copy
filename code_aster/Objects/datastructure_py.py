# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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


class OnlyParallelObject:
    """This object is only available in parallel."""

    def __init__(self, *args, **kwargs):
        raise NameError(
            "The object '{0}' is only available in parallel "
            "executions.".format(self.__class__.__name__)
        )


class PyDataStructure:
    """Temporary object used as a DataStructure during a Command execution."""

    def __init__(self, name="unnamed"):
        """Initialization"""
        self._name = name

    def getName(self):
        """Return the CO name."""
        return self._name

    @property
    def userName(self):
        """Same as 'getName'."""
        return self.getName()

    @userName.setter
    def userName(self, name):
        self._name = name

    def getType(self):
        """Return a type for syntax checking."""
        raise NotImplementedError("must be subclassed")


class AsInteger(PyDataStructure):
    """This class defines a simple integer used as a DataStructure."""

    @classmethod
    def getType(cls):
        """Return type as string."""
        return "ENTIER"


class AsFloat(PyDataStructure):
    """This class defines a simple float used as a DataStructure."""

    @classmethod
    def getType(cls):
        """Return type as string."""
        return "REEL"
