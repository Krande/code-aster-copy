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
Definition of options/flags for execution.
"""


class Options:
    """Enumerator for execution options.

    Some options are enabled and/or disabled from command line options.
    Some others are internally used by programming.

    - *Debug*: Debug mode.
    - *Abort*: Abort instead of raising an exception in case of errors.
    - *Continue*: Restart from an existing database.
    - *StrictUnpickling*: Fail when an object can not be unpickled.
    - *UseLegacyMode*: Create 'CO' objects instead of namedtuple.
    - *ShowDeprecated*: Show deprecation warnings.
    - *ShowSyntax*: Show syntax of commands.
    - *ShowChildCmd*: Show children commands called in depth.
    - *TestMode*: Testcase mode.
    - *SlaveMode*: Execution embedded by another program (do not abort, do not
      exit in case of time limit...).
    - *LastStep*: Last step of a study, database won't be reloaded from the
      temporary directory
    - *HPCMode*: High Performance Computing mode, parallel computing using
      domain decomposition.
    """

    Null = 0x0000
    Debug = 0x0001
    Abort = 0x0002
    Continue = 0x0004
    StrictUnpickling = 0x0008
    UseLegacyMode = 0x0010
    ShowDeprecated = 0x0020
    ShowSyntax = 0x0040
    ShowChildCmd = 0x0080
    TestMode = 0x0100
    SlaveMode = 0x0200
    LastStep = 0x0400
    HPCMode = 0x0800
    # do not forget to document each new option

    @classmethod
    def by_name(cls, name):
        """Return an option value by its name.
        *AttributeError* is raised if the option does not exist.

        Arguments:
            name (str): Option name.

        Returns:
            int: Option value.
        """
        return getattr(cls, name)
