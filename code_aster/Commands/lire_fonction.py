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

# person_in_charge: mathieu.courtois@edf.fr


from ..Messages import UTMESS
from ..Supervis import ExecuteMacro
from ..Utilities import deprecate


class ReadFunctionProperly(ExecuteMacro):
    """Command LIRE_FUNCTION that reads file values representing a function.
    """
    command_name = "LIRE_FONCTION"

    def adapt_syntax(self, keywords):
        """Hook to adapt syntax from a old version or for compatibility reasons.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords, changed
                in place.
        """
        # SEPAR is deprecated but tolerated syntax replaced by SEPARATEUR
        separ = keywords.pop("SEPAR", None)
        if separ:
            deprecate("LIRE_FONCTION/SEPAR", case=3,
                      help="Prefer use SEPARATEUR=... instead.")
            if "SEPARATEUR" not in keywords:
                if separ in ("None", ",", ";", "/"):
                    keywords["SEPARATEUR"] = separ
                else:
                    UTMESS('F', 'FONCT0_26')
            else:
                UTMESS('F', 'FONCT0_11')

LIRE_FONCTION = ReadFunctionProperly.run
