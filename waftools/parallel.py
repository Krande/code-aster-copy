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

import os
import os.path as osp
from functools import partial
from waflib import Configure, Errors, Utils, Logs


def options(self):
    self.load("compiler_c")
    self.load("compiler_cxx")
    self.load("compiler_fc")

    group = self.get_option_group("code_aster options")
    group.add_option(
        "--disable-mpi", dest="parallel", action="store_false", help="Build a sequential version"
    )
    group.add_option(
        "--enable-mpi",
        dest="parallel",
        action="store_true",
        default=os.environ.get("ENABLE_MPI") != "0",
        help="Build a parallel version with mpi (same as " "ENABLE_MPI environment variable)",
    )
    group.add_option(
        "--enable-openmp",
        dest="openmp",
        action="store_true",
        help="Build a parallel version supporting OpenMP",
    )
    group.add_option("--disable-openmp", dest="openmp", action="store_false", help="Disable OpenMP")
    group.add_option(
        "--enable-proc-status",
        dest="procstatus",
        action="store_true",
        help="force control of used memory with VmSize",
    )
    group.add_option(
        "--disable-proc-status",
        dest="procstatus",
        action="store_false",
        help="disable control of used memory with VmSize",
    )
    group.add_option(
        "--enable-ccache",
        dest="ccache",
        action="store_true",
        default=os.environ.get("ENABLE_CCACHE") == "1",
        help="Build using ccache if it is available, same as ENABLE_CCACHE=1",
    )
    group.add_option(
        "--disable-ccache",
        dest="ccache",
        action="store_false",
        help="Disable the use of ccache (this is the default)",
    )


def configure(self):
    opts = self.options
    # Configure.find_program uses first self.environ, then os.environ
    if opts.parallel:
        self.environ.setdefault("CC", "mpicc")
        self.environ.setdefault("CXX", "mpicxx")
        self.environ.setdefault("FC", "mpif90")
    else:
        self.environ.setdefault("CC", "gcc")
        self.environ.setdefault("CXX", "g++")
        self.environ.setdefault("FC", "gfortran")
    self.load_compilers()
    self.check_compilers_version()
    self.check_fortran_verbose_flag()
    self.check_openmp()
    # self.check_vmsize() is executed after mpiexec checking


###############################################################################


@Configure.conf
def check_compilers_version(self):
    self.start_msg("Checking for C compiler version")
    self.end_msg(self.env.CC_NAME.lower() + " " + ".".join(Utils.to_list(self.env.CC_VERSION)))
    # CXX_VERSION does not exist, c++ == c
    self.start_msg("Checking for Fortran compiler version")
    self.end_msg(self.env.FC_NAME.lower() + " " + ".".join(Utils.to_list(self.env.FC_VERSION)))


@Configure.conf
def load_compilers(self):
    if self.options.ccache:
        try:
            self.find_program("ccache")
            env = self.environ
            if "ccache" not in env["CC"]:
                env["CC"] = "ccache " + env["CC"]
            if "ccache" not in env["CXX"]:
                env["CXX"] = "ccache " + env["CXX"]
            if "ccache" not in env["FC"]:
                env["FC"] = "ccache " + env["FC"]
        except Errors.ConfigurationError:
            pass
    self.load("compiler_c")
    self.load("compiler_cxx")
    self.load("compiler_fc")
    if self.options.parallel:
        self.load_compilers_mpi()


@Configure.conf
def load_compilers_mpi(self):
    check = partial(
        self.check_cfg,
        args="--showme:compile --showme:link -show",
        package="",
        uselib_store="MPI",
        mandatory=False,
    )
    # raise ValueError(str(self.env.CC) + "  ///  " + str(self.env.CXX))

    ifort = "ifort" in self.env.FC_NAME.lower()
    icc = "icc" in self.env.CC_NAME.lower()

    # We won't alter environment if Intel compiler is detected...
    if not icc:
        msg = "Checking C compiler package (collect configuration flags)"
        if not check(path=self.env.CC, msg=msg):
            self.fatal("Unable to configure the parallel environment for C compiler")
        self.env["CCLINKFLAGS_MPI"] = self.env["LINKFLAGS_MPI"]
        del self.env["LINKFLAGS_MPI"]

    # We won't alter environment if Intel compiler is detected...
    if not ifort:
        msg = "Checking Fortran compiler package (collect configuration flags)"
        if not check(path=self.env.FC, msg=msg):
            self.fatal("Unable to configure the parallel environment for FORTRAN compiler")
        self.env["FCLINKFLAGS_MPI"] = self.env["LINKFLAGS_MPI"]
        del self.env["LINKFLAGS_MPI"]

    if not icc:
        self.check_cc(header_name="mpi.h", use="MPI")

    self.define("ASTER_HAVE_MPI", 1)
    self.env.BUILD_MPI = 1
    self.env.ASTER_HAVE_MPI = 1


@Configure.conf
def check_openmp(self):
    opts = self.options
    if opts.openmp is False:
        self.msg("Checking for OpenMP flag", "no", color="YELLOW")
        return
    # OpenMP interoperability is not secure
    # we consider both compiler should be from same vendor
    # Define CFLAGS_x and CCFLAGS_x to avoid ambiguous behaviour
    ifort = "ifort" in self.env.FC_NAME.lower()
    icc = "icc" in self.env.CC_NAME.lower()
    if ifort:
        self.env["FCFLAGS_OPENMP"] = ["/Qopenmp"]
        self.env["FCLINKFLAGS_OPENMP"] = ["/Qopenmp"]
        self.env["CFLAGS_OPENMP"] = self.env["FCFLAGS_OPENMP"]
        self.env["CCFLAGS_OPENMP"] = self.env["FCFLAGS_OPENMP"]
        self.env["CCLINKFLAGS_OPENMP"] = self.env["FCLINKFLAGS_OPENMP"]
        self.env["CXXFLAGS_OPENMP"] = self.env["FCFLAGS_OPENMP"]
        self.env["CXXLINKFLAGS_OPENMP"] = self.env["FCLINKFLAGS_OPENMP"]
        self.env.ASTER_HAVE_OPENMP = 1
        self.msg("Checking for OpenMP flag -qopenmp for Intel compilers", "yes", color="GREEN")
    elif not (ifort or icc):
        for x in ("-fopenmp", "-openmp", "-mp", "-xopenmp", "-omp", "-qsmp=omp"):
            try:
                self.check_fc(
                    msg="Checking for OpenMP flag %s" % x,
                    fragment="program main\n  call omp_get_num_threads()\nend program main",
                    fcflags=x,
                    fclinkflags=x,
                    uselib_store="OPENMP",
                )
            except self.errors.ConfigurationError:
                pass
            else:
                self.env["CFLAGS_OPENMP"] = self.env["FCFLAGS_OPENMP"]
                self.env["CCFLAGS_OPENMP"] = self.env["FCFLAGS_OPENMP"]
                self.env["CCLINKFLAGS_OPENMP"] = self.env["FCLINKFLAGS_OPENMP"]
                self.env["CXXFLAGS_OPENMP"] = self.env["FCFLAGS_OPENMP"]
                self.env["CXXLINKFLAGS_OPENMP"] = self.env["FCLINKFLAGS_OPENMP"]
                break
        else:
            self.fatal("Could not set OpenMP")
    else:
        self.fatal("Could not set OpenMP due to incompatible compilers...")
    self.define("ASTER_HAVE_OPENMP", 1)
    self.env.BUILD_OPENMP = 1


@Configure.conf
def check_sizeof_mpi_int(self):
    """Check size of MPI_Fint"""
    if self.get_define("ASTER_HAVE_MPI"):
        fragment = "\n".join(
            [
                "#include <stdio.h>",
                '#include "mpi.h"',
                "int main(void){",
                "    MPI_Fint var;",
                '    printf("%d", (int)sizeof(var));',
                "    return 0;",
                "}",
                "",
            ]
        )
        self.code_checker(
            "ASTER_MPI_INT_SIZE",
            self.check_cc,
            fragment,
            "Checking size of MPI_Fint integers",
            "unexpected value for sizeof(MPI_Fint): %(size)s",
            into=(4, 8),
            use="MPI",
        )


@Configure.conf
def check_vmsize(self):
    """Check for VmSize 'bug' with MPI or not."""
    opts = self.options
    if opts.procstatus is None:
        flag = self.get_define("ASTER_ENABLE_PROC_STATUS")
    else:
        flag = int(opts.procstatus)
    if flag is not None:
        self.start_msg("Check measure of VmSize using /proc")
        if flag not in (0, 1, "0", "1"):
            raise Errors.ConfigurationError("unexpected value: ASTER_ENABLE_PROC_STATUS=%s" % flag)
        flag = int(flag)
        if flag:
            self.end_msg("ok (ASTER_ENABLE_PROC_STATUS=%s)" % flag)
        else:
            self.end_msg("disabled (ASTER_ENABLE_PROC_STATUS=%s)" % flag, "YELLOW")
    elif not self.get_define("ASTER_HAVE_MPI"):
        self.start_msg("Check measure of VmSize using /proc")
        flag = 1
        self.end_msg("default (use /proc/PID/status)")
    else:
        self.start_msg("Checking measure of VmSize during MPI_Init")
        try:
            prg = osp.join(self.bldnode.abspath(), "test_mpi_init_" + str(os.getpid()))
            self.check_cc(fragment=fragment_failure_vmsize, mandatory=True, use="MPI", target=prg)
            cmd = self.env["base_mpiexec"] + ["-n", "1", prg]
            size = self.cmd_and_log(cmd)
        except Errors.WafError:
            self.end_msg(
                "failed (memory consumption can not be estimated " "during the calculation)",
                "YELLOW",
            )
        else:
            self.end_msg("ok (%s)" % size)
            flag = 1
            self.check_require_mpiexec(prg)
    if flag:
        self.define("ASTER_ENABLE_PROC_STATUS", 1)
    else:
        self.undefine("ASTER_ENABLE_PROC_STATUS")


@Configure.conf
def check_require_mpiexec(self, program):
    """Check if mpiexec is required to run a program with one process."""
    cfg = self.env["CONFIG_PARAMETERS"]
    previous = cfg["require_mpiexec"]
    cmt = ""
    self.start_msg("Checking if mpiexec is required")
    try:
        # run simple program without mpiexec
        self.cmd_and_log(program)
        required = 0
    except Errors.WafError:
        self.end_msg("mpiexec is required")
        required = 1
    else:
        if previous:
            cmt = ", but enabled by configuration parameters"
        self.end_msg("mpiexec is not required" + cmt)
    finally:
        os.remove(program)
    cfg["require_mpiexec"] = required or previous
    self.start_msg(". use 'require_mpiexec'")
    self.end_msg(str(cfg["require_mpiexec"]))


fragment_failure_vmsize = r"""
/*
   Check for unexpected value of VmPeak passing MPI_Init.
*/

#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

#include "mpi.h"


long read_vmsize()
{
    static char filename[80];
    static char sbuf[1024];
    char* str;
    int fd, size;

    sprintf(filename, "/proc/%ld/status", (long)getpid());
    fd = open(filename, O_RDONLY, 0);
    if (fd == -1)
        return -1;
    size = read(fd, sbuf, (sizeof sbuf) - 1);
    close(fd);

    str = strstr(sbuf, "VmSize:") + 8;
    return atol(str);
}

int main(int argc, char *argv[])
{
    int iret;
    long size;

    size = read_vmsize();

    MPI_Init(&argc, &argv);

    sleep(1);
    size = read_vmsize() - size;

    printf("%ld kB", size);
    if ( size > 10 * 1024 * 1024 ) {
        // it should be around 100 MB (dynlibs loaded)
        fprintf(stderr, "MPI initialization used more than 10 GB (%ld kB)\n", size);
        iret = 1;
    }
    else {
        iret = 0;
    }

    MPI_Finalize();

    return iret;
}
"""
