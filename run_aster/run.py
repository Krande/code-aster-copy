# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
:py:mod:`run` --- Main classes for execution
--------------------------------------------

This module defines the objects that prepare the working directory, copy the
data files, execute code_aster and copy the result files.
"""

import os
import os.path as osp
import tempfile
from glob import glob
from pathlib import Path
from subprocess import PIPE, run

from .command_files import add_import_commands, file_changed, stop_at_end
from .config import CFG
from .logger import WARNING, logger
from .status import StateOptions, Status, get_status
from .timer import Timer
from .utils import (
    RUNASTER_ROOT,
    RUNASTER_PLATFORM,
    cmd_abspath,
    compress,
    copy,
    make_writable,
    run_command,
    uncompress,
)

EXITCODE_FILE = "_exit_code_"
TMPMESS = "fort.6"
FMT_DIAG = """
------------------------------------------------------------
--- DIAGNOSTIC JOB : {state}
------------------------------------------------------------
"""


def create_temporary_dir(dir):
    """Create a temporary directory.

    Returns:
        str: Path of the directory.
    """
    os.makedirs(dir, exist_ok=True)
    return tempfile.mkdtemp(prefix="run_aster_", dir=dir)


class RunAster:
    """Execute code_aster from a `.export` file.

    Arguments:
        export (Export): Export object defining the calculation.
    """

    _show_comm = True

    @classmethod
    def factory(
        cls, export, test=False, env=False, tee=False, output=None, interactive=False, exectool=None
    ):
        """Return a *RunAster* object from an *Export* object.

        Arguments:
            export (Export): Export object defining the calculation.
            test (bool): for a testcase,
            env (bool): to only prepare the working directory and show
                command lines to be run,
            tee (bool): to follow execution output,
            interactive (bool): to keep Python interpreter active.
            exectool (str): command that preceeds code_aster command line.
            output (str): Path to redirect stdout.
        """
        class_ = RunAster
        if env:
            class_ = RunOnlyEnv
        return class_(export, test, tee, output, interactive, exectool)

    def __init__(
        self, export, test=False, tee=False, output=None, interactive=False, exectool=None
    ):
        self.export = export
        self.jobnum = str(os.getpid())
        logger.debug("Export content: %s", self.export.filename)
        logger.debug("\n%s", self.export)
        self._parallel = CFG.get("parallel", 0)
        self._test = test
        self._tee = tee
        self._output = output or TMPMESS
        self._interact = interactive
        self._exectool = exectool
        if self.export.get("hide-command"):
            self._show_comm = False
        procid = 0
        if self._parallel:
            procid = get_procid()
        procid = max(procid, 0)
        self._procid = procid
        # multi-steps study
        if not export.has_param("step"):
            export.set("step", 0)
        if not export.has_param("nbsteps"):
            export.set("nbsteps", len(self.export.commfiles))
        self._last = export.get("step") + 1 == export.get("nbsteps")

    def execute(self, wrkdir):
        """Execution in a working directory.

        Arguments:
            wrkdir (str): Working directory.

        Returns:
            Status: Status object.
        """
        if self._parallel:
            wrkdir = osp.join(wrkdir, f"proc.{self._procid}")
        os.makedirs(wrkdir, exist_ok=True)
        os.chdir(wrkdir)
        status = self._execute()
        return status

    def _execute(self):
        """Execution in the current working directory.

        Returns:
            Status: Status object.
        """
        timer = Timer()
        timer.load("__timer__")
        logger.info("TITLE Execution of code_aster")
        timer.start("Preparation of environment")
        self.prepare_current_directory()
        timer.stop()
        timer.start("Execution of code_aster")
        status = self.execute_study()
        timer.stop()
        if self._last or not status.is_completed():
            timer.start("Copying results")
            self.ending_execution(status.is_completed())
            logger.info("TITLE Execution summary")
            logger.info(timer.report())
            if self._procid == 0:
                logger.info(FMT_DIAG.format(state=status.diag))
        timer.save("__timer__")
        return status

    def prepare_current_directory(self):
        """Prepare the working directory."""
        logger.info("TITLE Prepare environment in %s", os.getcwd())
        self.export.write_to(self.jobnum + ".export")
        os.makedirs("REPE_IN", exist_ok=True)
        os.makedirs("REPE_OUT", exist_ok=True)

    def execute_study(self):
        """Execute the study.

        Returns:
            Status: Status object.
        """
        commfiles = [obj.path for obj in self.export.commfiles]
        if not commfiles:
            logger.error("no .comm file found")
        if self.export.get("nbsteps") > 1:
            os.makedirs("BASE_PREC", exist_ok=True)

        timeout = self.export.get("time_limit", 86400) * 1.25 * self.export.get("ncpus", 1)
        status = Status()
        comm = commfiles[0]
        logger.info(
            "TITLE Command file #{0} / {1}".format(
                self.export.get("step") + 1, self.export.get("nbsteps")
            )
        )
        comm = change_comm_file(comm, interact=self._interact, show=self._show_comm)
        status.update(self._exec_one(comm, timeout - status.times[-1]))
        self._coredump_analysis()
        return status

    def _exec_one(self, comm, timeout):
        """Show instructions for a command file.

        Arguments:
            comm (str): Command file name.
            timeout (float): Remaining time.

        Returns:
            Status: Status object.
        """
        idx = self.export.get("step")
        logger.info("TITLE Command line #%d:", idx + 1)
        timeout = int(max(1, timeout))
        cmd = self._get_cmdline(idx, comm, timeout)
        logger.info("    %s", " ".join(cmd))

        # use environment variable to make it works with ipython
        os.environ["PYTHONFAULTHANDLER"] = "1"
        exitcode = run_command(cmd, exitcode_file=EXITCODE_FILE)
        status = self._get_status(exitcode)
        msg = f"\nEXECUTION_CODE_ASTER_EXIT_{self.jobnum}={status.exitcode}\n\n"
        logger.info(msg)
        self._log_mess(msg)

        if status.is_completed():
            if not self._last:
                for vola in glob("vola.*"):
                    os.remove(vola)
                logger.info("saving result databases to 'BASE_PREC'...")
                for base in glob("glob.*") + glob("bhdf.*") + glob("pick.*"):
                    copy(base, "BASE_PREC")
            msg = f"execution ended (command file #{idx + 1}): {status.diag}"
            logger.info(msg)
        else:
            logger.info("restoring result databases from 'BASE_PREC'...")
            for base in glob(osp.join("BASE_PREC", "*")):
                copy(base, os.getcwd())
            msg = f"execution failed (command file #{idx + 1}): {status.diag}"
            logger.warning(msg)
        if self._procid == 0:
            self._log_mess(FMT_DIAG.format(state=status.diag))
        return status

    def _get_cmdline_exec(self, commfile, idx):
        """Build the command line really executed, without redirection.

        Arguments:
            commfile (str): Command file name.
            idx (int): Index of execution.

        Returns:
            list[str]: List of command line arguments, without redirection.
        """
        cmd = []
        if self._exectool:
            cmd.append(self._exectool)
        wrapped = False
        python = CFG.get("python")
        if not self._interact:
            # absolute path is necessary to call the debugger
            cmd.append(cmd_abspath(python))
            if self._parallel:
                # see documentation of `mpi4py.run`
                cmd.extend(["-m", "mpi4py"])
        else:
            cmd.append(CFG.get("python_interactive", python))
            wrapped = CFG.get("python_interactive_is_wrapped")
            print("DEBUG:", cmd[-1], wrapped)
            cmd.append("-i")
        # To show executed lines with trace module:
        # import sys
        # ign = [sys.prefix, sys.exec_prefix, "$HOME/.local", os.getenv("PYTHONPATH")]
        # cmd.extend(["-m", "trace", "--trace",
        #             "--ignore-dir=" + ":".join(ign)])
        cmd.append(commfile)
        if wrapped:
            cmd.append("--")
        # remaining arguments are treated by code_aster script
        if self._test:
            cmd.append("--test")
        if self._last:
            cmd.append("--last")
        # copy datafiles only the first time because all share the same workdir
        if idx == 0:
            for obj in self.export.datafiles:
                cmd.append(f'--link="{obj.as_argument}"')
        cmd.extend(self.export.args)
        # TODO add pid + mode to identify the process by asrun
        return cmd

    def _coredump_analysis(self):
        """Process the coredump file."""
        core = glob("core*")
        if not core:
            return
        logger.info("\ncoredump analysis...")
        python3 = cmd_abspath(CFG.get("python"))
        if not osp.isfile(python3):
            logger.warning("'python3' not found in PATH.")
            return
        tmpf = "cmd_gdb.sh"
        with open(tmpf, "w") as fobj:
            fobj.write(os.linesep.join(["where", "quit", ""]))
        # may be required if gdb is linked against a different libpython
        bck = os.environ.pop("PYTHONHOME", None)
        cmd = ["gdb", "-batch", "-x", tmpf, "-e", python3, "-c", core[0]]
        run_command(cmd)
        if bck:
            os.environ["PYTHONHOME"] = bck

    def _get_cmdline(self, idx, commfile, timeout):
        """Build the command line.

        Arguments:
            idx (int): Index of execution.
            commfile (str): Command file name.
            timeout (float): Remaining time.

        Returns:
            list[str]: List of command line arguments.
        """
        cmd = self._get_cmdline_exec(commfile, idx)
        if self._interact:
            if not Path(self._output).exists():
                Path(self._output).touch()
        elif self._tee:
            orig = " ".join(cmd)
            if RUNASTER_PLATFORM == "linux":
                cmd = [f"( {orig} ; echo $? > {EXITCODE_FILE} )", "2>&1"]
            else:
                cmd = [f"( {orig} & echo %errorlevel% > {EXITCODE_FILE} )"]
            cmd.extend(["|", "tee", "-a", self._output])
        else:
            cmd.extend([">>", self._output, "2>&1"])
        if RUNASTER_PLATFORM == "linux":
            cmd.insert(0, f"ulimit -c unlimited ; ulimit -t {timeout:.0f} ;")
        return cmd

    def _get_status(self, exitcode):
        """Get the execution status.

        Arguments:
            exitcode (int): Return code.

        Returns:
            Status: Status object.
        """
        status = get_status(exitcode, self._output, test=self._test and self._last)
        expected = self.export.get("expected_diag", [])
        if status.diag in expected:
            status.state = StateOptions.Ok
            status.exitcode = 0
        elif expected:
            status.update(Status(StateOptions.Fatal, 1))
        return status

    def ending_execution(self, is_completed):
        """Post execution phase : copying results, cleanup...

        Arguments:
            is_completed (bool): *True* if execution succeeded,
                *False* otherwise.
        """
        logger.info("TITLE Content of %s after execution:", os.getcwd())
        logger.info(_ls(".", "REPE_OUT"))
        if self._procid != 0:
            return

        results = self.export.resultfiles
        if results:
            logger.info("TITLE Copying results")
            copy_resultfiles(results, is_completed, test=self._test)

    def _log_mess(self, msg):
        """Log a message into the *message* file."""
        with open(self._output, "a") as fobj:
            fobj.write(msg + "\n")


class RunOnlyEnv(RunAster):
    """Prepare a working directory for a manual execution.

    Arguments:
        export (Export): Export object defining the calculation.
    """

    _show_comm = False

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self._procid > 0:
            logger.setLevel(WARNING)

    def prepare_current_directory(self):
        """Prepare the working directory."""
        if self.export.get("step") == 0:
            super().prepare_current_directory()

    def execute_study(self):
        """Execute the study.

        Returns:
            Status: Status object.
        """
        if self.export.get("step") == 0:
            logger.info("TITLE Copy/paste these command lines:")
            profile = osp.join(RUNASTER_ROOT, "share", "aster", "profile.sh")
            timeout = self.export.get("time_limit", 0) * 1.25
            if not self._parallel:
                logger.info("    cd %s", os.getcwd())
            else:
                logger.info("    cd %s", osp.dirname(os.getcwd()))
            logger.info("    . %s", profile)
            logger.info("    ulimit -c unlimited")
            logger.info("    ulimit -t %.0f", timeout)
        return super().execute_study()

    def _exec_one(self, comm, timeout):
        """Show instructions for a command file.

        Arguments:
            comm (str): Command file name.
            timeout (float): Remaining time.
        """
        idx = self.export.get("step")
        cmd = []
        if self._parallel:
            cmd.append("#!/bin/bash")
            cmd.append("cd proc.$({0})".format(CFG.get("mpi_get_rank")))
        cmd.append(" ".join(self._get_cmdline_exec(comm, idx)))
        cmd.append("")
        shell = f"cmd{idx}.sh"
        with open(shell, "w") as fobj:
            fobj.write("\n".join(cmd))
        os.chmod(f"cmd{idx}.sh", 0o755)
        if not self._parallel:
            logger.info(cmd[0])
        elif self._procid == 0:
            args_cmd = dict(mpi_nbcpu=self.export.get("mpi_nbcpu", 1), program="proc.0/" + shell)
            cmd = CFG.get("mpiexec").format(**args_cmd)
            logger.info(cmd)
        return Status(StateOptions.Ok, exitcode=0)

    def ending_execution(self, _):
        """Nothing to do in this case."""


def get_procid():
    """Return the identifier of the current process.

    Returns:
        int: Process ID, -1 if a parallel version is not run under *mpiexec*.
    """
    proc = run(CFG.get("mpi_get_rank"), shell=True, stdout=PIPE, universal_newlines=True)
    try:
        procid = int(proc.stdout.strip())
    except ValueError:
        procid = -1
    return procid


def get_nbcores():
    """Return the number of available cores.

    Returns:
        int: Number of cores.
    """
    return os.cpu_count()


def change_comm_file(comm, interact=False, wrkdir=None, show=False):
    """Change a command file.

    Arguments:
        comm (str): Command file name.
        wrkdir (str, optional): Working directory to write the changed file
            if necessary (defaults: current working directory).
        show (bool): Show file content if *True*.

    Returns:
        str: Name of the file to be executed (== *comm* if nothing changed)
    """
    with open(comm, "rb") as fobj:
        text_init = fobj.read().decode(errors="replace")
    text = add_import_commands(text_init)
    if interact:
        text = stop_at_end(text)
    changed = text.strip() != text_init.strip()
    if changed:
        text = file_changed(text, comm)
    if show:
        logger.info("\nContent of the file to execute:\n%s\n", text)
    if not changed:
        return comm

    filename = osp.join(wrkdir or ".", osp.basename(comm) + ".changed.py")
    with open(filename, "wb") as fobj:
        fobj.write(text.encode())
    return filename


def copy_datafiles(files, verbose=True):
    """Copy data files into the working directory.

    Arguments:
        files (list[File]): List of File objects.
        verbose (bool): Verbosity.
    """
    for obj in files:
        dest = None
        # fort.*
        if obj.unit != 0 or obj.filetype == "nom":
            if obj.unit in (6, 15):
                raise IndexError(
                    "Files fort.6 and fort.15 are reserved.\n"
                    "Please change unit number for: {}".format(obj.path)
                )
            dest = "fort." + str(obj.unit)
            if obj.filetype == "nom":
                dest = osp.basename(obj.path)
            # warning if file already exists
            if osp.exists(dest):
                logger.warning("%r overwrites %r", obj.path, dest)
            if obj.compr:
                dest += ".gz"
        # for directories
        else:
            if obj.filetype in ("base", "bhdf"):
                dest = osp.basename(obj.path)
            elif obj.filetype == "repe":
                dest = "REPE_IN"

        if dest is not None:
            copy(obj.path, dest, verbose=verbose)
            if obj.compr:
                dest = uncompress(dest)
            # move the bases in main directory
            if obj.filetype in ("base", "bhdf"):
                for fname in glob(osp.join(dest, "*")):
                    os.rename(fname, osp.basename(fname))
            # force the file to be writable
            make_writable(dest)


def copy_resultfiles(files, copybase, test=False):
    """Copy result files from the working directory.

    Arguments:
        files (list[File]): List of File objects.
        copybase (bool): Tell if result databases will be copied.
        test (bool, optional): *True* for a testcase, *False* for a study.
    """
    for obj in files:
        if test and obj.unit not in (6, 15):
            continue
        lsrc = []
        # fort.*
        if obj.unit != 0:
            lsrc.append("fort." + str(obj.unit))
        elif obj.filetype == "nom":
            lsrc.append(osp.basename(obj.path))
        # for directories
        else:
            if copybase and obj.filetype in ("base", "bhdf"):
                lbase = glob("bhdf.*")
                if not lbase or obj.filetype == "base":
                    lbase = glob("glob.*")
                lsrc.extend(lbase)
                lsrc.extend(glob("pick.*"))
            elif obj.filetype == "repe":
                lsrc.extend(glob(osp.join("REPE_OUT", "*")))

        for filename in lsrc:
            if not osp.exists(filename):
                logger.warning("file not found: %s", filename)
            else:
                if obj.compr:
                    filename = compress(filename)
                if obj.isdir and not osp.exists(obj.path):
                    os.makedirs(obj.path)
                copy(filename, obj.path, verbose=True)


def _ls(*paths):
    if RUNASTER_PLATFORM == "linux":
        proc = run(["ls", "-l"] + list(paths), stdout=PIPE, universal_newlines=True)
    else:
        proc = run(["dir"] + list(paths), stdout=PIPE, universal_newlines=True, shell=True)
    return proc.stdout
