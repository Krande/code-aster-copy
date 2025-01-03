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

# person_in_charge: mathieu.courtois at edf.fr

"""Ce module définit la classe `SYSTEM` et la fonction `ExecCommand`
qui est présente uniquement pour commodité pour les Macros.

La classe SYSTEM est semblable à celle utilisée dans asrun.
"""

__all__ = ["SYSTEM", "ExecCommand"]

import io
import os
import re
import sys
import tempfile

from .i18n import translate as _


def _exitcode(status, default=0):
    """
    Extract the exit code from status. Return 'default' if the process
    has not been terminated by exit.
    """
    if os.WIFEXITED(status):
        iret = os.WEXITSTATUS(status)
    elif os.WIFSIGNALED(status):
        iret = os.WTERMSIG(status)
    elif os.WIFSTOPPED(status):
        iret = os.WSTOPSIG(status)
    else:
        iret = default
    return iret


class SYSTEM:

    """
    Class to encapsultate "system" commands (this a simplified version of
    ASTER_SYSTEM class defined in ASTK_SERV part).
    """

    # this value should be set during installation step.
    MaxCmdLen = 1024
    # line length -9
    _LineLen = 80 - 9

    def __init__(self, **kargs):
        """
        Initialization.
        Optionnal arguments : silent, verbose, debug, cc_files, maxcmdlen.
        """
        self.verbose = kargs.get("verbose", True)
        self.debug = kargs.get("debug", False)
        self.cc_files = kargs.get("cc_files", None)
        if "maxcmdlen" in kargs:
            self.MaxCmdLen = kargs["maxcmdlen"]

    def _mess(self, msg, cod=""):
        """
        Just print a message
        """
        self._print("%-18s %s" % (cod, msg))

    def _print(self, *args, **kargs):
        """
        print replacement.
        Optionnal argument :
           term  : line terminator (default to os.linesep).
        """
        term = kargs.get("term", os.linesep)
        files = set([sys.stdout])
        if self.cc_files:
            files.add(self.cc_files)
        for f in files:
            if isinstance(f, io.IOBase):
                txt = " ".join(["%s" % a for a in args])
                f.write(txt.replace(os.linesep + " ", os.linesep) + term)
                f.flush()
            else:
                print(_("file object expected : %s / %s") % (type(f), repr(f)))

    def VerbStart(self, cmd, verbose=None):
        """
        Start message in verbose mode
        """
        Lm = self._LineLen
        if verbose is None:
            verbose = self.verbose
        if verbose:
            pcmd = cmd
            if len(cmd) > Lm - 2 or cmd.count("\n") > 0:
                pcmd = pcmd + "\n" + " " * Lm
            self._print(("%-" + str(Lm) + "s") % (pcmd,), term="")

    def VerbEnd(self, iret, output="", verbose=None):
        """
        End message in verbose mode
        """
        if verbose is None:
            verbose = self.verbose
        if verbose:
            if iret == 0:
                self._print("[  OK  ]")
            else:
                self._print(_("[FAILED]"))
                self._print(_("Exit code : %d") % iret)
            if (iret != 0 or self.debug) and output:
                self._print(output)

    def VerbIgnore(self, verbose=None):
        """
        End message in verbose mode
        """
        if verbose is None:
            verbose = self.verbose
        if verbose:
            self._print(_("[ SKIP ]"))

    def Shell(
        self,
        cmd,
        bg=False,
        verbose=False,
        follow_output=False,
        alt_comment=None,
        interact=False,
        separated_stderr=False,
    ):
        """
        Execute a command shell
           cmd           : command
           bg            : put command in background if True
           verbose       : print status messages during execution if True
           follow_output : follow interactively output of command
           alt_comment   : print this "alternative comment" instead of "cmd"
           interact      : allow the user to interact with the process.
        Return :
           iret     : exit code if bg = False,
                      0 if bg = True
           output   : output lines (as string)
        """
        if not alt_comment:
            alt_comment = cmd
        if bg:
            interact = False
        if len(cmd) > self.MaxCmdLen:
            self._mess(
                (_("length of command shell greater " "than %d characters.") % self.MaxCmdLen),
                "<A>_ALARM",
            )
        if self.debug:
            self._print("cmd :", cmd, "background : %s" % bg, "follow_output : %s" % follow_output)
        self.VerbStart(alt_comment, verbose=verbose)
        if follow_output and verbose:
            print(os.linesep + _("Command output :"))

        fout = tempfile.NamedTemporaryFile(delete=False)
        fout.close()
        ferr = tempfile.NamedTemporaryFile(delete=False)
        ferr.close()
        if bg:
            new_cmd = cmd + " &"
        elif follow_output:
            new_cmd = "( %s ) | tee %s" % (cmd, fout.name)
        else:
            if not separated_stderr:
                new_cmd = "( %s ) > %s 2>&1" % (cmd, fout.name)
            else:
                new_cmd = "( %s ) > %s 2> %s" % (cmd, fout.name, ferr.name)
        if self.debug:
            self._print("modified cmd :", new_cmd)
        # execution
        iret = os.system(new_cmd)
        fout = open(fout.name, "rb")
        output = fout.read().decode()
        ferr = open(ferr.name, "rb")
        error = ferr.read().decode()
        fout.close()
        ferr.close()

        if follow_output:
            # repeat header message
            self.VerbStart(alt_comment, verbose=verbose)
        mat = re.search("EXIT_CODE=([0-9]+)", output)
        if mat:
            iret = int(mat.group(1))
        self.VerbEnd(iret, output, verbose=verbose)
        if self.debug and iret != 0:
            self._print("<debug> ERROR : iret = %s" % iret)
            self._print("STDOUT", output, all=True)
            self._print("STDERR", error, all=True)
        if bg:
            iret = 0
        if not separated_stderr:
            result = iret, output
        else:
            result = iret, output, error
        return result


# Juste par commodité.
system = SYSTEM()
ExecCommand = system.Shell


if __name__ == "__main__":
    iret, output = ExecCommand("ls", alt_comment="Lancement de la commande...")
