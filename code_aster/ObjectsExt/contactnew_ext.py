# coding=utf-8
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

"""
:py:class:`ContactNew` --- Assignment of contact
************************************************************************
"""

from libaster import ContactNew

from ..Utilities import injector

@injector(ContactNew)
class ExtendedContactNew:
    cata_sdj = "SD.sd_char_contact.sd_char_contact"

    def __getinitargs__(self):
        """Returns the argument required to reinitialize a
        ContactNew object during unpickling.
        """
        return (self.getName(), self.getModel())
