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

import os
import os.path as osp
import re
from functools import partial

from waflib import Configure, Errors, Utils


def options(self):
    self.load("compiler_c")
    self.load("compiler_cxx")
    self.load("compiler_fc")

    group = self.get_option_group("code_aster options")
    group.add_option(
        "--disable-mpi",
        dest="parallel",
        default=os.environ.get("ENABLE_MPI") != "0",
        action="store_false",
        help="Build a sequential version",
    )
    group.add_option(
        "--enable-mpi",
        dest="parallel",
        action="store_true",
        help="Build a parallel version with mpi (same as ENABLE_MPI environment variable)",
    )
    group.add_option(
        "--enable-openmp",
        dest="openmp",
        action="store_true",
        default=None,
        help="Build a version supporting OpenMP",
    )
    group.add_option("--disable-openmp", dest="openmp", action="store_false", help="Disable OpenMP")
    group.add_option(
        "--use-srun",
        dest="use_srun",
        default=os.environ.get("ENABLE_SRUN") == "1",
        action="store_true",
        help="use 'srun' to start processes instead of 'mpiexec' (default: False), "
        "can be set with environment variable ENABLE_SRUN=1",
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
    self.load_compilers()
    self.check_fortran_verbose_flag()
    self.check_openmp()
    # self.check_vmsize() is executed after mpiexec checking


###############################################################################


@Configure.conf
def check_compilers_type_version(self):
    def _get_version(prog):
        output = self.cmd_and_log(Utils.to_list(prog) + ["--version"])
        return output.splitlines()[0]

    self.env.CC_IS_INTEL = False
    self.env.FC_IS_INTEL = False
    self.start_msg("Checking for C compiler version")
    line1 = _get_version(self.env.CC)
    expr = re.compile("(intel|oneapi|icx|icc)")
    if expr.search(line1.lower()):
        self.env.CC_IS_INTEL = True
    self.end_msg(line1)
    # CXX_VERSION does not exist, c++ == c
    self.start_msg("Checking for Fortran compiler version")
    line1 = _get_version(self.env.FC)
    expr = re.compile("(intel|oneapi|ifx|ifort)")
    if expr.search(line1.lower()):
        self.env.FC_IS_INTEL = True
    self.end_msg(line1)


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
    self.check_compilers_type_version()
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

    # We won't alter environment if Intel compiler is detected...
    if not self.env.CC_IS_INTEL:
        msg = "Checking C compiler package (collect configuration flags)"
        if not check(path=self.env.CC, msg=msg):
            self.fatal("Unable to configure the parallel environment for C compiler")
        self.env["CCLINKFLAGS_MPI"] = self.env["LINKFLAGS_MPI"]
        del self.env["LINKFLAGS_MPI"]

    # We won't alter environment if Intel compiler is detected...
    if not self.env.FC_IS_INTEL:
        msg = "Checking Fortran compiler package (collect configuration flags)"
        if not check(path=self.env.FC, msg=msg):
            self.fatal("Unable to configure the parallel environment for FORTRAN compiler")
        self.env["FCLINKFLAGS_MPI"] = self.env["LINKFLAGS_MPI"]
        del self.env["LINKFLAGS_MPI"]

    self.define("ASTER_HAVE_MPI", 1)
    self.env.BUILD_MPI = 1
    self.env.ASTER_HAVE_MPI = 1
    self.check_mpi_fortran_interface()


@Configure.conf
def check_openmp(self):
    opts = self.options
    if opts.openmp is False:
        self.msg("Checking for OpenMP flag", "disabled", color="YELLOW")
        return
    # OpenMP interoperability is not secure
    # we consider both compiler should be from same vendor
    # Define CFLAGS_x and CCFLAGS_x to avoid ambiguous behaviour
    if self.env.FC_IS_INTEL and self.env.CC_IS_INTEL:
        self.env["FCFLAGS_OPENMP"] = ["-qopenmp"]
        self.env["FCLINKFLAGS_OPENMP"] = ["-qopenmp"]
        self.env["CFLAGS_OPENMP"] = self.env["FCFLAGS_OPENMP"]
        self.env["CCFLAGS_OPENMP"] = self.env["FCFLAGS_OPENMP"]
        self.env["CCLINKFLAGS_OPENMP"] = self.env["FCLINKFLAGS_OPENMP"]
        self.env["CXXFLAGS_OPENMP"] = self.env["FCFLAGS_OPENMP"]
        self.env["CXXLINKFLAGS_OPENMP"] = self.env["FCLINKFLAGS_OPENMP"]
        self.env.ASTER_HAVE_OPENMP = 1
        self.msg("Checking for OpenMP flag -qopenmp for Intel compilers", "yes", color="GREEN")
    elif not (self.env.FC_IS_INTEL or self.env.CC_IS_INTEL):
        for x in ("-fopenmp", "-qopenmp", "-openmp", "-mp", "-xopenmp", "-omp", "-qsmp=omp"):
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
    if not self.get_define("ASTER_HAVE_MPI"):
        return
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
def check_mpi_fortran_interface(self):
    """Check mpi fortran module"""
    self.check_cc(header_name="mpi.h", use="MPI")
    # because use of mpif.h is deprecated (does not work with flang)
    fragment = "\n".join(["program main", "use mpi_f08", "end program main", ""])
    self.check_fc(
        fragment=fragment,
        msg="Checking for mpi_f08 module",
        errmsg="no, but not yet used!",
        use="MPI",
        mandatory=False,
    )


@Configure.conf
def check_vmsize(self):
    """Check for VmSize 'bug' with MPI or not."""
    if not self.get_define("ASTER_HAVE_MPI"):
        return
    self.start_msg("Checking measure of VmSize during MPI_Init")
    try:
        prg = osp.join(self.bldnode.abspath(), "test_mpi_init_" + str(os.getpid()))
        self.check_cc(fragment=fragment_failure_vmsize, mandatory=True, use="MPI", target=prg)
        cmd = self.env["base_mpiexec"] + ["-n", "1", prg]
        size = self.cmd_and_log(cmd)
    except Errors.WafError:
        raise Errors.ConfigurationError(
            "failed (memory consumption can not be estimated during the calculation)"
        )
    self.end_msg("ok (%s)" % size)
    self.check_require_mpiexec(prg)


@Configure.conf
def check_require_mpiexec(self, program):
    """Check if mpiexec/srun is required to run a program with one process."""
    cfg = self.env["CONFIG_PARAMETERS"]
    previous = cfg["require_mpiexec"]
    starter = "mpiexec" if "srun" not in program else "srun"
    cmt = ""
    self.start_msg(f"Checking if {starter} is required")
    try:
        # run simple program without mpiexec
        self.cmd_and_log(program)
        required = 0
    except Errors.WafError:
        self.end_msg(f"{starter} is required")
        required = 1
    else:
        if previous:
            cmt = ", but enabled by configuration parameters"
        self.end_msg(f"{starter} is not required" + cmt)
    finally:
        os.remove(program)
    cfg["require_mpiexec"] = required or previous
    self.start_msg(". use 'require_mpiexec'")
    self.end_msg(str(cfg["require_mpiexec"]))


fragment_failure_vmsize = r"""
/*
   Check for unexpected value of VmPeak passing MPI_Init.
   Extracted from mempid.c
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
    pid_t pid = getpid();

    sprintf(filename, "/proc/%ld/status", pid);
    fd = open(filename, O_RDONLY, 0);
    if (fd == -1)
        return -1;
    size = read(fd, sbuf, (sizeof sbuf) - 1);
    close(fd);

    str = strstr(sbuf, "VmSize:") + 8;
    return atol(str);
}

long read_memory_currrent() {
    char cgroup_path[256] = { 0 };
    char memory_path[512] = { 0 };
    FILE *cgroup_file = NULL;
    FILE *memory_file = NULL;
    long memory_bytes;
    pid_t pid = getpid();

    snprintf( cgroup_path, sizeof( cgroup_path ), "/proc/%d/cgroup", pid );
    cgroup_file = fopen( cgroup_path, "r" );
    if ( !cgroup_file ) {
        return -1;
    }

    if ( !fgets( cgroup_path, sizeof( cgroup_path ), cgroup_file ) ) {
        fclose( cgroup_file );
        return -1;
    }
    fclose( cgroup_file );

    cgroup_path[strcspn( cgroup_path, "\n" )] = '\0';
    char *cgroup_rel_path = strchr( cgroup_path, '/' );
    if ( !cgroup_rel_path ) {
        return -1;
    }

    // memory.current
    snprintf( memory_path, sizeof( memory_path ), "/sys/fs/cgroup%s/memory.current",
              cgroup_rel_path );
    memory_file = fopen( memory_path, "r" );
    if ( !memory_file ) {
        return -1;
    }
    if ( fscanf( memory_file, "%ld", &memory_bytes ) != 1 ) {
        fclose( memory_file );
        return -1;
    }
    fclose( memory_file );
    return memory_bytes / 1024;
}

#if defined ASTER_ENABLE_CGROUP_STATS
#define get_vmsize read_memory_currrent
#elif defined ASTER_ENABLE_PROC_STATUS
#define get_vmsize read_vmsize
#else
#error ERROR: needs ASTER_ENABLE_PROC_STATUS or ASTER_ENABLE_CGROUP_STATS
#endif

int main(int argc, char *argv[])
{
    int iret;
    long size;

    size = get_vmsize();

    MPI_Init(&argc, &argv);

    sleep(1);
    size = get_vmsize() - size;

    printf("%ld kB", size);
#if defined ASTER_ENABLE_CGROUP_STATS
    printf(", using 'read_memory_currrent'");
#elif defined ASTER_ENABLE_PROC_STATUS
    printf(", using 'read_vmsize'");
#endif
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
