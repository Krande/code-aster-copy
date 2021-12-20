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
:py:class:`DEBUT` --- Initialization of code_aster
**************************************************

The :py:class:`Starter` starts the execution by initializing the code_aster
memory manager (*Jeveux*). For this task, it parses the arguments through the
:py:class:`~code_aster.Utilities.ExecutionParameter.ExecutionParameter` object.
By default, arguments are those read from the command line and those passed
:py:func:`.init`.
If the command line arguments are not set for code_aster, they can be ignored
with ``code_aster.init(..., noargv=True)``.

Some Python objects that have to be available from :py:mod:`libaster` are
passed during the initialization to the
:py:class:`~code_aster.Utilities.ExecutionParameter.ExecutionParameter`.
"""

import aster_core
import libaster
from run_aster.run import copy_datafiles

from ..Behaviours import catalc
from ..Cata.Syntax import tr
from ..Cata.SyntaxUtils import remove_none
from ..Helpers import LogicalUnitFile
from ..Messages import UTMESS, MessageLog
from ..Supervis import CommandSyntax, ExecuteCommand, Serializer, loadObjects
from ..Supervis.code_file import track_coverage
from ..Supervis.ctopy import checksd, print_header
from ..Supervis.TestResult import testresu_print
from ..Utilities import MPI, ExecutionParameter, Options, import_object, logger
from ..Utilities.i18n import localization

try:
    import debugpy

    HAS_DEBUGPY = True
except ImportError:
    HAS_DEBUGPY = False


class ExecutionStarter:
    """Initialize the
    :class:`~code_aster.Utilities.ExecutionParameter.ExecutionParameter` object
    for requests from the both sides Python/Fortran."""

    params = _is_initialized = None

    @classmethod
    def init(cls, argv=None, fcomm=0):
        """Initialization of class attributes.

        Attributes:
            argv (list[str]): List of command line arguments.
            fcomm (int, optional): Id of the MPI communicator.

        Returns:
            bool: *True* if the initialization has been done, *False* if the
            execution was already initialized.
        """
        if cls._is_initialized:
            return False
        params = cls.params = ExecutionParameter()
        params.parse_args(argv)
        params.catalc = catalc
        params.logical_unit = LogicalUnitFile
        params.syntax = CommandSyntax
        params.print_header = print_header
        params.checksd = checksd
        params.testresu_print = testresu_print
        copy_datafiles(params.export.datafiles)
        aster_core.register(params, MessageLog)
        libaster.jeveux_init(fcomm)
        cls._is_initialized = True
        return True


class Starter(ExecuteCommand):
    """Define the command DEBUT."""

    command_name = "DEBUT"
    arg_init = []

    @staticmethod
    def _code_enabled(keywords):
        """Tell if CODE is enabled.

        Arguments:
            keywords (dict): User's keywords, changed in place.

        Returns:
            bool: *True* if CODE is present, *False* otherwise.
        """
        return keywords.get("CODE")

    @classmethod
    def run(cls, **keywords):
        """Run the Command.

        Arguments:
            keywords (dict): User keywords
        """
        if not ExecutionStarter.init(cls.arg_init):
            return

        super(Starter, cls).run(**keywords)

    @classmethod
    def _run_with_argv(cls, **keywords):
        """Wrapper to have the same depth calling loadObjects..."""
        cmd = cls()
        cmd._result = None
        cmd._cata.addDefaultKeywords(keywords)
        remove_none(keywords)
        cmd.exec_(keywords)

    @classmethod
    def run_with_argv(cls, **keywords):
        """Run the command with the arguments from the command line.

        Arguments:
            keywords (dict): User keywords
        """
        cls._run_with_argv(**keywords)

    def exec_(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        iwarn = False
        stop_with = "EXCEPTION"
        if ExecutionParameter().option & Options.Abort:
            stop_with = "ABORT"
        if ExecutionParameter().option & Options.TestMode or self._code_enabled(
            keywords
        ):
            ExecutionParameter().enable(Options.TestMode)
            stop_with = "ABORT"
            iwarn = True
            track_coverage(self._cata, self.command_name, keywords)

        erreur = keywords.get("ERREUR")
        if erreur:
            if erreur.get("ERREUR_F"):
                stop_with = erreur["ERREUR_F"]
        if ExecutionParameter().option & Options.SlaveMode:
            stop_with = "EXCEPTION"
        # must be the first call to correctly set 'vini' in onerrf
        libaster.onFatalError(stop_with)

        debug = keywords.get("DEBUG")
        if debug:
            jxveri = debug.get("JXVERI", "NON") == "OUI"
            ExecutionParameter().set_option("jxveri", int(jxveri))
            if jxveri:
                UTMESS("I", "SUPERVIS_23")
            sdveri = debug.get("SDVERI", "NON") == "OUI"
            ExecutionParameter().set_option("sdveri", int(sdveri))
            if sdveri:
                UTMESS("I", "SUPERVIS_24")
            dbgjeveux = (
                debug.get("JEVEUX", "NON") == "OUI"
                or debug.get("VERI_BASE") is not None
            )
            ExecutionParameter().set_option("dbgjeveux", int(dbgjeveux))
            if dbgjeveux:
                UTMESS("I", "SUPERVIS_12")
            iwarn = iwarn or jxveri or sdveri or dbgjeveux
        if iwarn:
            UTMESS("I", "SUPERVIS_22", valk=("--test", "code_aster.init()"))
        if ExecutionParameter().get_option("hook_post_exec"):
            path = ExecutionParameter().get_option("hook_post_exec")
            hook = import_object(path)
            self.register_hook(hook)

        if keywords.get("IMPR_MACRO") == "OUI":
            ExecutionParameter().enable(Options.ShowChildCmd)

        if keywords.get("LANG"):
            translation = localization.translation(keywords["LANG"])
            tr.set_translator(translation.gettext)

        if keywords.get("IGNORE_ALARM"):
            for idmess in keywords["IGNORE_ALARM"]:
                MessageLog.disable_alarm(idmess)
        super(Starter, self).exec_(keywords)

    def _call_oper(self, syntax):
        """Call fortran operator.

        Arguments:
            syntax (*CommandSyntax*): Syntax description with user keywords.
        """
        logger.info("starting the execution...")
        libaster.call_debut(syntax)


class Restarter(Starter):
    """Define the command POURSUITE."""

    command_name = "POURSUITE"
    arg_init = ["--continue"]

    @staticmethod
    def _code_enabled(keywords):
        """Tell if CODE is enabled.

        Arguments:
            keywords (dict): User's keywords, changed in place.

        Returns:
            bool: *True* if CODE is present, *False* otherwise.
        """
        return keywords.get("CODE") == "OUI"

    def _call_oper(self, syntax):
        """Call fortran operator.

        Arguments:
            syntax (*CommandSyntax*): Syntax description with user keywords.
        """
        if not Serializer.canRestart():
            logger.error("restart aborted!")

        logger.info("restarting from a previous execution...")
        libaster.call_poursuite(syntax)
        # 1:_call_oper, 2:ExecuteCommand.exec_, 3:Starter.exec_,
        #  4:Restarter.run, 5:ExecuteCommand.run_, 6:ExecuteCmd.run, 7:user
        # 1:_call_oper, 2:ExecuteCommand.exec_, 3:Starter.exec_,
        #  4:_run_with_argv, 5:run_with_argv, 6:init, 7:user
        loadObjects(level=7)


DEBUT = Starter.run
POURSUITE = Restarter.run


def init(*argv, **kwargs):
    """Initialize code_aster as `DEBUT`/`POURSUITE` command does + command
    line options.

    If the code_aster study is embedded under another Python program, the
    "--slave" option may be useful to catch exceptions (even *TimeLimitError*)
    and **try** not to exit the Python interpreter.

    Arguments:
        argv (list): List of command line arguments.
        kwargs (dict): Keywords arguments passed to 'DEBUT'/'POURSUITE' +
            'debug (bool)' same as -g/--debug,
            'debugpy (int)' to start a debugpy session on this port number,
            'noargv (bool)' to ignore previously passed arguments,
            'comm (*mpi4py.MPI.Comm*)' to select the MPI communicator.
    """
    if kwargs.pop("noargv", False):
        ExecutionParameter().set_argv([])
    comm = kwargs.pop("comm", None)
    fcomm = comm.py2f() if comm else 0
    if not ExecutionStarter.init(argv, fcomm):
        return

    if kwargs.pop("debug", False):
        ExecutionParameter().enable(Options.Debug)

    if kwargs.get("debugpy") and HAS_DEBUGPY:
        debugpy.listen(("localhost", kwargs["debugpy"]))
        print("Waiting for debugger attach")
        debugpy.wait_for_client()
        debugpy.breakpoint()
        # add 10 hours for debugging
        tpmax = ExecutionParameter().get_option("tpmax")
        ExecutionParameter().set_option("tpmax", tpmax + 36000)
    kwargs.pop("debugpy", None)

    if ExecutionStarter.params.option & Options.Continue:
        Restarter.run_with_argv(**kwargs)
    else:
        Starter.run_with_argv(**kwargs)
