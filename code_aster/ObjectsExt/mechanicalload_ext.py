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

# person_in_charge: mathieu.courtois@edf.fr
"""
:py:class:`MechanicalLoad` --- Assignment of mechanical load
*******************************************************************
"""

from libaster import MechanicalLoadComplex, MechanicalLoadFunction, MechanicalLoadReal

from ..Utilities import injector


@injector(MechanicalLoadReal)
class ExtendedMechanicalLoadReal:
    cata_sdj = "SD.sd_char_meca.sd_char_meca"

    def __getinitargs__(self):
        """Returns the argument required to reinitialize a MechanicalLoadReal
        object during unpickling.
        """
        return (self.getName(), self.getModel())


@injector(MechanicalLoadFunction)
class ExtendedMechanicalLoadFunction:
    cata_sdj = "SD.sd_char_meca.sd_char_meca"

    def __getinitargs__(self):
        """Returns the argument required to reinitialize a MechanicalLoadFunction
        object during unpickling.
        """
        return (self.getName(), self.getModel())


@injector(MechanicalLoadComplex)
class ExtendedMechanicalLoadComplex:
    cata_sdj = "SD.sd_char_meca.sd_char_meca"

    def __getinitargs__(self):
        """Returns the argument required to reinitialize a MechanicalLoadComplex
        object during unpickling.
        """
        return (self.getName(), self.getModel())
