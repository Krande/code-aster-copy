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

from ..Objects import GenericMechanicalLoad
from ..Supervis import ExecuteCommand
from ..Utilities import deprecate, force_list

class MechanicalLoadDefinition(ExecuteCommand):
    """Command that creates the
    :class:`~code_aster.Objects.GenericMechanicalLoad`"""
    command_name = "AFFE_CHAR_MECA_C"
    def adapt_syntax(self, keywords):
        """Adapt keywords.
        Replace LIAISON = 'ENCASTRE'           
        """
        common_dofs = ('DX', 'DY', 'DZ', 'DRX', 'DRY', 'DRZ')
        # replace DDL_IMPO/LIAISON by DDL_IMPO/DX=0
        keywords["DDL_IMPO"] = force_list(keywords.get("DDL_IMPO", []))
        for fact in keywords["DDL_IMPO"]:
            block = fact.pop("LIAISON", None)
            if block == 'ENCASTRE':
                deprecate("DLL_IMPO/LIAISON='ENCASTRE'", case=3, level=5,
                          help="Use BLOCAGE = ('DEPLACEMENT','ROTATION')")
                for ddl in common_dofs:
                    fact[ddl] = 0.
    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        self._result = GenericMechanicalLoad(keywords["MODELE"])


AFFE_CHAR_MECA_C = MechanicalLoadDefinition.run
