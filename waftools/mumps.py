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
import os
import pathlib
import re
from functools import partial

from waflib import Configure, Errors, Utils, Logs


def options(self):
    group = self.add_option_group("Mumps library options")
    group.add_option(
        "--disable-mumps",
        action="store_false",
        default=None,
        dest="enable_mumps",
        help="Disable MUMPS support",
    )
    group.add_option(
        "--enable-mumps",
        action="store_true",
        default=None,
        dest="enable_mumps",
        help="Force MUMPS support",
    )
    group.add_option(
        "--mumps-libs",
        type="string",
        dest="mumps_libs",
        default=None,
        help="mumps librairies to use when linking",
    )
    group.add_option(
        "--embed-mumps",
        dest="embed_mumps",
        default=False,
        action="store_true",
        help="Embed MUMPS libraries as static library",
    )


def configure(self):
    if self.env.CC_NAME == "msvc":
        mumps_seq_incl = pathlib.Path(os.environ.get("CONDA_PREFIX")) / 'Library/include/mumps_seq'
        self.env.append_value('INCLUDES', [mumps_seq_incl.as_posix()])

    try:
        self.env.stash()
        self.check_mumps()
    except Errors.ConfigurationError as exc:
        self.reset_msg()
        self.env.revert()
        self.undefine("ASTER_HAVE_MUMPS")
        if self.options.enable_mumps is True:
            raise Errors.ConfigurationError(
                str(exc) + "\n" + "check your environment or use --disable-mumps"
            )
    else:
        self.define("ASTER_HAVE_MUMPS", 1)


###############################################################################
@Configure.conf
def check_mumps(self):
    opts = self.options
    if opts.enable_mumps is False:
        raise Errors.ConfigurationError("MUMPS disabled")
    self.check_mumps_headers()
    try:
        self.check_mumps_version(("5.6", "5.5", "5.4"))
    except Errors.ConfigurationError:
        self.check_mumps_version(("5.7",), beta=True)
    self.check_sizeof_mumps_integer()
    if opts.mumps_libs is None:
        opts.mumps_libs = "dmumps zmumps smumps cmumps mumps_common pord"
    if not self.env.BUILD_MPI and self.env.CC_NAME != "msvc":
        opts.mumps_libs += " mpiseq"
    if opts.mumps_libs:
        self.check_mumps_libs()


@Configure.conf
def check_mumps_libs(self):
    opts = self.options
    check_mumps = partial(
        self.check_fc, uselib_store="MUMPS", use="MUMPS MPI MATH OPENMP", mandatory=True
    )
    if opts.embed_all or opts.embed_mumps:
        check = lambda lib: check_mumps(stlib=lib)
    else:
        check = lambda lib: check_mumps(lib=lib)
    list(map(check, Utils.to_list(opts.mumps_libs)))


@Configure.conf
def check_mumps_headers(self):
    fragment = r"""
      PROGRAM MAIN
      INCLUDE '{0}'
      PRINT *, 'ok'
      END PROGRAM MAIN
"""

    # Logs.warn(f"{self.env['INCLUDES']=}")
    headers = [i + "mumps_struc.h" for i in "sdcz"] + ["mpif.h"]
    if self.get_define("ASTER_HAVE_MPI"):
        for path in self.env["INCLUDES"][:]:
            if "include_seq" in path:
                self.env["INCLUDES"].remove(path)
                self.start_msg("Removing path from INCLUDES")
                self.end_msg(path, "YELLOW")
    for inc in headers:
        try:

            self.start_msg("Checking for {0}".format(inc))
            self.check_fc(
                fragment=fragment.format(inc),
                compile_filename="test.F",
                uselib_store="MUMPS",
                uselib="MUMPS MPI",
                mandatory=True,
            )
        except:
            self.end_msg("no", "YELLOW")
            raise
        else:
            self.end_msg("yes")


@Configure.conf
def check_mumps_version(self, expected_versions, beta=False):
    fragment = r"""
#include <stdio.h>
#include "smumps_c.h"

int main(void){
    printf("%s", MUMPS_VERSION);
    return 0;
}"""
    self.start_msg("Checking mumps version")
    vers = ""
    try:
        ret = self.check_cc(
            fragment=fragment, use="MUMPS", mandatory=True, execute=True, define_ret=True
        )
        self.env["MUMPS_VERSION"] = ret
        vers = re.sub("(consortium|c| ) *$", "", ret)
        if vers[:3] not in expected_versions:
            raise Errors.ConfigurationError("expected versions: {0}".format(expected_versions))

        self.define("ASTER_MUMPS_REDUCMPI", 1)
        if re.search("(consortium|c) *$", ret):
            self.define("ASTER_MUMPS_CONSORTIUM", 1)
    except:
        if vers:
            vers = " (%s)" % self.env["MUMPS_VERSION"]
        self.end_msg("no" + vers, "YELLOW")
        raise
    else:
        self.define("ASTER_MUMPS_VERSION", ret)
        self.end_msg(self.env["MUMPS_VERSION"] + (" beta support" if beta else ""))


@Configure.conf
def check_sizeof_mumps_integer(self):
    """Check Mumps integer size
    include "dmumps_struc.h"
    type(dmumps_struc) :: dmpsk
    print *, sizeof(i)
    end
    """
    fragment = "\n".join(
        [
            'include "dmumps_struc.h"',
            "type(dmumps_struc) :: dmpsk",
            "print *, sizeof(dmpsk%n)",
            "end",
        ]
    )
    self.code_checker(
        "ASTER_MUMPS_INT_SIZE",
        self.check_fc,
        fragment,
        "Checking size of Mumps integer",
        "unexpected value for sizeof(mumps_int): %(size)s",
        into=(4, 8),
        use="MUMPS",
    )
