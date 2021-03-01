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
import gc

from ..Objects import DataStructure, Function
from ..Utilities import deprecate, force_list, get_caller_context
from ..Supervis import ExecuteCommand


class Deleter(ExecuteCommand):
    """Command that deletes *DataStructure* instances from the calling stack."""
    command_name = "DETRUIRE"

    def compat_syntax(self, keywords):
        """Hook to adapt syntax from a old version or for compatibility reasons.

        Arguments:
            keywords (dict): User's keywords, changed in place.
        """
        to_del = []
        kwlist = force_list(keywords.pop("CONCEPT", []))
        if kwlist:
            deprecate("DETRUIRE/CONCEPT/NOM", case=4, level=5,
                      help="Just use DETRUIRE/NOM=... instead.")
            for occ in kwlist:
                to_del.extend(force_list(occ["NOM"]))
            keywords["NOM"] = to_del
        if keywords.pop("OBJET", None):
            deprecate("DETRUIRE/OBJET", case=4, level=5,
                      help="Use DETRUIRE/NOM=... instead.")
            # to have at least one object
            keywords.setdefault("NOM", Function())

    def exec_(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords, changed in place to force
                deletion. "NOM" will not be available for 'post_exec'.
        """
        if self.level > 1:
            deprecate("DETRUIRE should be used in a macro-command", case=9,
                      level=5)
            return

        to_del = [obj.getName() for obj in keywords.pop("NOM")]
        if not to_del:
            return

        # calling context
        context = get_caller_context(3)

        for name in list(context.keys()):
            if isinstance(context[name], DataStructure):
                if context[name].getName() in to_del:
                    del context[name]

        # force garbage collection
        gc.collect()


DETRUIRE = Deleter.run
