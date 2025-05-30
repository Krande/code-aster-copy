# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
import re
from functools import partial
from waflib import Configure, Utils, Errors


def options(self):
    group = self.add_option_group("Metis library options")
    group.add_option(
        "--disable-metis",
        action="store_false",
        default=None,
        dest="enable_metis",
        help="Disable METIS support",
    )
    group.add_option(
        "--enable-metis",
        action="store_true",
        default=None,
        dest="enable_metis",
        help="Force METIS support",
    )
    group.add_option(
        "--metis-libs",
        type=str,
        dest="metis_libs",
        default=None,
        help="metis librairies to use when linking",
    )
    group.add_option(
        "--embed-metis",
        dest="embed_metis",
        default=False,
        action="store_true",
        help="Embed METIS libraries as static library",
    )


def configure(self):
    try:
        self.env.stash()
        self.check_metis()
    except Errors.ConfigurationError:
        self.reset_msg()
        self.env.revert()
        self.undefine("ASTER_HAVE_METIS")
        if self.options.enable_metis:
            raise
    else:
        self.define("ASTER_HAVE_METIS", 1)


###############################################################################
@Configure.conf
def check_metis(self):
    opts = self.options
    if opts.enable_metis is False:
        raise Errors.ConfigurationError("METIS disabled")
    self.check_metis_libs()
    self.check_metis_headers()
    self.check_metis_idx()
    self.check_metis_real()
    self.check_metis_version()


@Configure.conf
def check_metis_libs(self):
    opts = self.options
    check_metis = partial(
        self.check_cc, uselib_store="METIS", use="PARMETIS METIS M", mandatory=True
    )
    if opts.embed_all or opts.embed_metis:
        check = lambda lib: check_metis(stlib=lib)
    else:
        check = lambda lib: check_metis(lib=lib)
    # METIS is currently provided by ParMETIS
    if opts.metis_libs is None:
        if not check_metis(lib="metis", mandatory=False):
            check_metis(lib="parmetis")
    else:
        list(map(check, Utils.to_list(opts.metis_libs)))


@Configure.conf
def check_metis_headers(self):
    if self.is_defined("ASTER_PLATFORM_MINGW"):
        self.define("USE_GKREGEX", 1)
    check = partial(
        self.check_cc, header_name="metis.h", uselib_store="METIS", use="PARMETIS METIS M"
    )
    self.start_msg("Checking for header metis.h")
    try:
        if not check(mandatory=False):
            if not check(includes=[osp.join(self.env.INCLUDEDIR, "metis")], mandatory=False):
                check(includes=[osp.join(self.env.OLDINCLUDEDIR, "metis")], mandatory=True)
    except:
        self.end_msg("no", "YELLOW")
        raise
    else:
        self.end_msg("yes")


@Configure.conf
def check_metis_idx(self):
    fragment = r"""
#include <stdio.h>
#include <metis.h>
int main(void){
    idx_t idx;
    printf("%d", (int)sizeof(idx));
    return 0;
}"""
    self.code_checker(
        "ASTER_METIS_IDX_SIZE",
        self.check_cc,
        fragment,
        "Checking size of metis idx_t",
        "unexpected value for sizeof(idx_t): %(size)s",
        into=(4, 8),
        use="PARMETIS METIS M",
    )


@Configure.conf
def check_metis_real(self):
    fragment = r"""
#include <stdio.h>
#include <metis.h>
int main(void){
    real_t real;
    printf("%d", (int)(sizeof(real)));
    return 0;
}"""
    self.code_checker(
        "ASTER_METIS_REAL_SIZE",
        self.check_cc,
        fragment,
        "Checking size of metis real_t",
        "unexpected value for sizeof(real_t): %(size)s",
        into=(4, 8),
        use="PARMETIS METIS M",
    )


@Configure.conf
def check_metis_version(self):
    fragment = r"""
#include <stdio.h>
#include <metis.h>
int main(void){
#ifdef METISTITLE
/* metis 4 */
    printf("METISTITLE: %s", METISTITLE);
    return 0;
#endif
#if defined(METIS_VER_MAJOR) && defined(METIS_VER_MINOR) && defined(METIS_VER_SUBMINOR)
    printf("METISVER: %d.%d.%d", METIS_VER_MAJOR, METIS_VER_MINOR, METIS_VER_SUBMINOR);
    return 0;
#endif
/* unexpected */
    return 1;
}"""
    self.start_msg("Checking metis version")
    try:
        ret = self.check_cc(
            fragment=fragment, use="PARMETIS METIS M", mandatory=True, execute=True, define_ret=True
        )
        mat4 = re.search(r"METISTITLE: *METIS *(?P<vers>[0-9]+\.[0-9]+\.\w+) ", ret)
        mat5 = re.search(r"METISVER: *(?P<vers>[0-9]+\.[0-9]+\.\w+)", ret)
        vers = (mat4 and mat4.group("vers")) or (mat5 and mat5.group("vers"))
        major = int(vers.split(".")[0])
        if major != 5:
            self.end_msg("unsupported metis version: %s" % vers, "RED")
            raise Errors.ConfigurationError
    except:
        self.end_msg("can not get version", "RED")
        raise
    else:
        self.end_msg(vers)
