# coding: utf-8

# Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

from ..Messages import UTMESS
from ..Objects import ElementaryCharacteristics
from ..Supervis import ExecuteCommand


class EltCharacteristicsAssignment(ExecuteCommand):

    """Command that assigns
    :class:`~code_aster.Objects.ElementaryCharacteristics` on a model.
    """

    command_name = "AFFE_CARA_ELEM"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        self._result = ElementaryCharacteristics(keywords["MODELE"])

    def adapt_syntax(self, keywords):
        """Hook to adapt syntax *after* syntax checking.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords, changed
                in place.
        """
        if "MASSIF" in keywords:
            # Check that MASSIF appears once only if there is TOUT in MASSIF simple keywords
            l_dic_kws = keywords.get("MASSIF")
            if type(l_dic_kws) == tuple:  # il y a plus d'une occurrence de MASSIF
                for dic in l_dic_kws:
                    if "TOUT" in dic.keys():
                        UTMESS("F", "SUPERVIS_10")


AFFE_CARA_ELEM = EltCharacteristicsAssignment.run
