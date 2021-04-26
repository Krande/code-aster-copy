# coding: utf-8

# Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

from ..Objects import Mesh
from ..Supervis import ExecuteCommand
from ..Utilities import deprecate, force_list


class MeshModification(ExecuteCommand):
    """Command that changes a :class:`~code_aster.Objects.Mesh`.
    """
    command_name = "MODI_MAILLAGE"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        self._result = keywords["MAILLAGE"]

    def add_dependencies(self, keywords):
        """Register input *DataStructure* objects as dependencies.

        Arguments:
            keywords (dict): User's keywords.
        """

    @staticmethod
    def compat_syntax(keywords):
        """ Update Keyword ORIE_PEAU_2D(3D) to ORIE_PEAU"""
        keywords_mapping = {
            'ORIE_PEAU_2D': {'GROUP_MA': 'GROUP_MA_PEAU', 'GROUP_MA_SURF': 'GROUP_MA_INTERNE'},
            'ORIE_PEAU_3D': {'GROUP_MA': 'GROUP_MA_PEAU', 'GROUP_MA_VOLU': 'GROUP_MA_INTERNE'}
        }
        for kw in keywords_mapping:
            kws = force_list(keywords.get(kw, []))
            if kws:
                deprecate("MODI_MAILLAGE/ORIE_PEAU_2D(3D)=_F(GROUP_MA=(...), GROUP_MA_SURF(VOLU)=(...))", case=3,
                          level=5,
                          help="Use MODI_MAILLAGE/ORIE_PEAU=_F(GROUP_MA_PEAU=(...), GROUP_MA_INTERNE=(...))")
                keywords['ORIE_PEAU'] = keywords.get(kw)
                keywords.pop(kw)
                for fact in kws:
                    keys = list(fact.keys())
                    for key in keys:
                        fact[keywords_mapping[kw][key]] = fact[key]
                        fact.pop(key)


MODI_MAILLAGE = MeshModification.run
