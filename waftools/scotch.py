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

import os.path as osp
from functools import partial
from waflib import Configure, Utils, Errors


def options(self):
    group = self.add_option_group("Scotch library options")
    group.add_option(
        "--disable-scotch",
        action="store_false",
        default=None,
        dest="enable_scotch",
        help="Disable SCOTCH support",
    )
    group.add_option(
        "--enable-scotch",
        action="store_true",
        default=None,
        dest="enable_scotch",
        help="Force scotch support",
    )
    group.add_option(
        "--scotch-libs",
        type="string",
        dest="scotch_libs",
        default=None,
        help="scotch librairies used when linking",
    )
    group.add_option(
        "--embed-scotch",
        dest="embed_scotch",
        default=False,
        action="store_true",
        help="Embed SCOTCH libraries as static library",
    )


def configure(self):
    try:
        self.env.stash()
        self.check_scotch()
    except Errors.ConfigurationError:
        self.reset_msg()
        self.env.revert()
        self.undefine("ASTER_HAVE_SCOTCH")
        self.undefine("ASTER_HAVE_PTSCOTCH")
        if self.options.enable_scotch is True:
            raise
    else:
        self.define("ASTER_HAVE_SCOTCH", 1)
        if self.env.BUILD_MPI:
            self.define("ASTER_HAVE_PTSCOTCH", 1)


###############################################################################


@Configure.conf
def check_scotch(self):
    opts = self.options
    if opts.enable_scotch is False:
        raise Errors.ConfigurationError("SCOTCH disabled")

    self.check_scotch_headers()
    self.check_scotch_version()

    if opts.scotch_libs is None:
        opts.scotch_libs = "esmumps scotch scotcherr"
        if self.env.BUILD_MPI:
            opts.scotch_libs += " ptscotch ptscotcherr"

    # code_aster v11.0.1: FICHE 016627
    if "scotchmetis" in opts.scotch_libs:
        raise Errors.ConfigurationError(
            "scotchmetis variant library is not " "compatible with code_aster"
        )

    if opts.scotch_libs:
        self.check_scotch_libs()


@Configure.conf
def check_scotch_libs(self):
    opts = self.options
    check_scotch = partial(
        self.check_cc, mandatory=True, uselib_store="SCOTCH", use="MPI", uselib="SCOTCH M Z"
    )
    if opts.embed_all or opts.embed_scotch:
        check_lib = lambda lib: check_scotch(stlib=lib)
    else:
        check_lib = lambda lib: check_scotch(lib=lib)
    check_lib(Utils.to_list(opts.scotch_libs))


@Configure.conf
def check_scotch_headers(self):
    self.start_msg("Checking for header scotch.h")
    headers = "stdio.h stdlib.h stdint.h sys/types.h scotch.h"
    # can not test 'ptscotch.h' without including 'mpi.h' before...
    # if self.env.BUILD_MPI:
    #     headers += " ptscotch.h"
    try:
        check = partial(
            self.check, header_name=headers, uselib_store="SCOTCH", use="MPI", uselib="SCOTCH MPI"
        )

        if not check(mandatory=False):
            if not check(includes=[osp.join(self.env.INCLUDEDIR, "scotch")], mandatory=False):
                check(includes=[osp.join(self.env.OLDINCLUDEDIR, "scotch")], mandatory=True)
    except:
        self.end_msg("no", "YELLOW")
        raise
    else:
        self.end_msg("yes")


@Configure.conf
def check_scotch_version(self):
    # scotch.h may use int64_t without including <sys/types.h>
    fragment = r"""
#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include "scotch.h"

int main(void){
    printf("(%d, %d, %d)", (int)SCOTCH_VERSION, (int)SCOTCH_RELEASE,
           (int)SCOTCH_PATCHLEVEL);
    return 0;
}"""
    self.start_msg("Checking scotch version")
    mvers = self.env["MUMPS_VERSION"] or "0"
    mumps_vers = tuple([int(num) for num in mvers.split(".")[:2]])
    required = 7 if mumps_vers >= (5, 5) else 6
    try:
        ret = self.check_cc(
            fragment=fragment,
            use="SCOTCH",
            uselib_store="SCOTCH",
            mandatory=True,
            execute=True,
            define_ret=True,
        )
        vers = eval(ret)
        if vers[0] < required:
            self.end_msg(
                "Scotch version >= {1} is required, unexpected version: {0}".format(ret, required),
                "RED",
            )
            raise Errors.ConfigurationError
    except:
        self.end_msg("no", "YELLOW")
        raise
    else:
        self.end_msg(".".join([str(i) for i in eval(ret)]))
