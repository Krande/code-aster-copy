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

from ..Cata.Syntax import _F
from ..Messages import UTMESS
from ..Supervis import ExecuteCommand
from ..Utilities import force_list


class ImprResu(ExecuteCommand):
    """Command IMPR_RESU.
    """
    command_name = "IMPR_RESU"

    def add_result_name(self, resu):
        """Try to add NOM_RESU_MED keyword if not set.

        Arguments:
            resu (dict): factor keyword occurrence of RESU, changed in place.
        """
        # if resu.get("NOM_RESU_MED") or resu.get("NOM_CHAM_MED"):
        #     return
        if resu.get("RESULTAT"):
            resu_name = resu["RESULTAT"].userName[0:8]
            if not resu_name:
                UTMESS('A', 'MED3_2', valk=resu["RESULTAT"].getName())
            elif resu.get("NOM_CHAM"):
                resu.pop("NOM_RESU_MED", None)
                resu_name = resu_name.ljust(8, "_")
                if not resu.get("NOM_CHAM_MED"):
                    resu["NOM_CHAM"] = force_list(resu["NOM_CHAM"])
                    resu["NOM_CHAM_MED"] = [resu_name + field
                                            for field in resu["NOM_CHAM"]]
            elif not resu.get("NOM_RESU_MED"):
                resu["NOM_RESU_MED"] = resu_name
        elif resu.get("CHAM_GD"):
            field_name = resu["CHAM_GD"].userName[0:8]
            if resu.get("NOM_CHAM_MED"):
                return
            if not field_name:
                UTMESS('A', 'MED3_2', valk=resu["CHAM_GD"].getName())
            else:
                resu["NOM_CHAM_MED"] = field_name


    def adapt_syntax(self, keywords):
        """Hook to adapt syntax from a old version or for compatibility reasons.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords, changed
                in place.
        """
        if keywords.get("FORMAT") in (None, "MED"):
            keywords["RESU"] = force_list(keywords["RESU"])
            for resu in keywords["RESU"]:
                self.add_result_name(resu)


IMPR_RESU = ImprResu.run
