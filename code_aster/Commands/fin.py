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

# person_in_charge: mathieu.courtois@edf.fr
"""
:py:class:`FIN` --- Finalization of code_aster
**********************************************

The :py:class:`Closer` finalizes the execution by closing the code_aster
memory manager (*Jeveux*).

The objects existing in the context where :py:func:`~code_aster.Commands.FIN`
is called are pickled while their *Jeveux* content is saved in ``glob.*``
databases.
A call to the function :py:func:`~code_aster.Supervis.Serializer.saveObjects`
is equivalent.
"""

import sys
import libaster

from ..Supervis import ExecuteCommand, saveObjects, FinalizeOptions as Options
from ..Utilities import ExecutionParameter, logger


class Closer(ExecuteCommand):
    """Command that closes the execution."""
    command_name = "FIN"
    _options = None
    _exit = None

    def compat_syntax(self, keywords):
        """Consume argument to force to exit after the command.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords, changed
                in place.
        """
        self._exit = keywords.pop("exit", False)
        # removed keyword
        keywords.pop("FORMAT_HDF", None)

    def exec_(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        self._options = Options.SaveBase
        if keywords.get("INFO_RESU") == "OUI":
            self._options |= Options.InfoResu
        if keywords.get("RETASSAGE") == "OUI":
            self._options |= Options.Repack
        if keywords.get("PROC0") == "OUI":
            self._options |= Options.OnlyProc0
        super().exec_(keywords)

    def _call_oper(self, dummy):
        """Save objects that exist in the context of the caller.

        The memory manager is closed when :mod:`libaster` is unloaded.
        """
        # Ensure that `saveObjects` has not been already called
        if libaster.jeveux_status():
            # 1: here, 2: exec_, 3: parent.exec_, 4: run_, 5: run, 6: user space
            saveObjects(level=6, options=self._options)

        logger.info(repr(ExecutionParameter().timer))

    def post_exec(self, keywords):
        """Force to exit if `exit=True` option was passed.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        if self._exit:
            sys.exit()


FIN = Closer.run
