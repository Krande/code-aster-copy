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
from itertools import product
from functools import partial
from subprocess import Popen, PIPE
from waflib import Configure, Errors, Logs, Utils

BLAS = ("flexiblas", "openblas", "blas")
BLACS = ("blacs",)
LAPACK = ("lapack",)
SCALAPACK = ("scalapack",)
OPTIONAL_DEPS = ("cblas",)


def options(self):
    group = self.add_option_group("Mathematics libraries options")
    group.add_option(
        "--maths-libs",
        type=str,
        dest="maths_libs",
        default=None,
        help="List of libraries that provide blas and lapack symbols. "
        'Otherwise, one may use special values: "auto" to search them automatically, '
        '"mkl", "flexiblas" or "openblas" to force a type of libraries.',
    )
    group.add_option(
        "--embed-maths",
        dest="embed_math",
        default=False,
        action="store_true",
        help="Embed math libraries as static library",
    )


def configure(self):
    # always check for libpthread, libm (never in static)
    self.check_cc(uselib_store="M", lib="m")
    self.check_cc(uselib_store="Z", lib="z")
    self.check_number_cores()
    if self.options.maths_libs == "mkl":
        if not self.detect_mkl():
            raise Errors.ConfigurationError(
                "can not find MKL libraries, try '--maths-libs=auto' instead."
            )
    elif self.options.maths_libs in BLAS:
        self.detect_math_lib([self.options.maths_libs])
    elif self.options.maths_libs in (None, "auto"):
        # try MKL first, then automatic blas/lapack
        if not self.detect_mkl():
            self.detect_math_lib()
    elif self.options.maths_libs:
        self.check_opts_math_lib()
        self.check_math_libs_call_blas_lapack("RED")
    self.check_libm_after_files()


###############################################################################
@Configure.conf
def check_opts_math_lib(self):
    opts = self.options
    embed = opts.embed_math or opts.embed_all

    def check_lib(lib):
        return self.check_cc(
            **{
                "mandatory": True,
                "uselib_store": "MATH",
                "use": "MATH MPI",
                ("st" * embed + "lib"): lib,
            }
        )

    for lib in Utils.to_list(opts.maths_libs):
        check_lib(lib)


@Configure.conf
def check_sizeof_blas_int(self):
    """Check size of blas integers"""
    self.set_define_from_env(
        "ASTER_BLAS_INT_SIZE",
        "Setting size of blas/lapack integers",
        "unexpected value for blas int: %(size)s",
        into=(4, 8),
        default=4,
    )


@Configure.conf
def check_libm_after_files(self):
    """Avoid warning #10315: specifying -lm before files may supercede the
    Intel(R) math library and affect performance"""
    self.start_msg("Setting libm after files")
    flags = self.env.LINKFLAGS_CLIB
    if "-lm" in flags:
        while True:
            try:
                flags.remove("-lm")
            except ValueError:
                break
        self.end_msg('yes ("-lm" removed from LINKFLAGS_CLIB)')
        self.env.LINKFLAGS_CLIB = flags
    else:
        self.end_msg("nothing done")


@Configure.conf
def detect_mkl(self):
    """Try to detect MKL"""
    opts = self.options
    # MKL can be installed either as a standalone package or with Intel
    # compiler. MKLROOT may be undefined (conda package for example)
    self.start_msg("Detecting MKL libraries")
    suffix = "_lp64" if "64" in self.env.DEST_CPU else ""
    scalapack = ""
    blacs = []
    thread = "mkl_sequential"
    core = "mkl_core"
    libs = []
    # http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/
    if "ifort" in self.env.FC_NAME.lower() or "icc" in self.env.CC_NAME.lower():
        if self.get_define("ASTER_HAVE_OPENMP"):
            thread = "mkl_intel_thread"
        interf = "mkl_intel" + suffix
        if self.get_define("ASTER_HAVE_MPI") and opts.enable_mumps:
            # scalapack = "mkl_scalapack" + suffix
            # blacs = "mkl_blacs_intelmpi" + suffix
            scalapack = "scalapack"
    else:
        if self.get_define("ASTER_HAVE_OPENMP"):
            thread = "mkl_gnu_thread"
        interf = "mkl_gf" + suffix
        if self.get_define("ASTER_HAVE_MPI") and opts.enable_mumps:
            # This needs to add all libs into LD_PRELOAD (libmpi.so + all mkl libs...)
            # scalapack = "mkl_scalapack" + suffix
            # blacs = "mkl_blacs_openmpi" + suffix
            scalapack = "scalapack"
    libs.append(interf)
    libs.append(thread)
    libs.append(core)
    if scalapack:
        libs.append(scalapack)
    if blacs:
        libs.append(blacs)
    self.end_msg("trying " + str(libs))
    try:
        self.env.stash()
        self.env.append_value("LIB_MATH", libs)
        if "MKLROOT" in os.environ:
            self.env.append_value("LIBPATH_MATH", os.environ["MKLROOT"] + "/lib/intel64")
        self.check_math_libs_call(color="YELLOW")
    except:
        self.env.revert()
        return False
    else:
        self.define("ASTER_HAVE_MKL", 1)
        self.env.commit()
        return True


@Configure.conf
def detect_math_lib(self, libs=BLAS):
    opts = self.options
    embed = opts.embed_math or (opts.embed_all and not self.get_define("ASTER_HAVE_MPI"))
    varlib = ("ST" if embed else "") + "LIB_MATH"

    # blas
    blaslibs, lapacklibs = self.get_mathlib_from_numpy()
    self.check_math_libs(list(libs) + blaslibs, embed)

    # lapack
    opt_lapack = False
    if "openblas" in self.env.get_flat(varlib) or "flexiblas" in self.env.get_flat(varlib):
        # check that lapack is embedded in openblas/flexiblas
        try:
            self.check_math_libs_call_blas_lapack(color="YELLOW")
            opt_lapack = True
        except:
            pass
    if not opt_lapack:
        self.check_math_libs(list(LAPACK) + lapacklibs, embed)
        self.check_math_libs_call_blas_lapack()

    def _scalapack():
        """Check scalapack"""
        libs = list(SCALAPACK)
        libs = libs + ["".join(n) for n in product(libs, ["mpi", "-mpi", "openmpi", "-openmpi"])]
        return self.check_math_libs(libs, embed)

    def _blacs():
        """Check blacs"""
        libs = list(BLACS)
        libs = (
            libs
            + ["".join(n) for n in product(libs, ["mpi", "-mpi", "openmpi", "-openmpi"])]
            + ["".join(n) for n in product(["mpi", "mpi-", "openmpi", "openmpi-"], libs)]
        )  # check the 3 blacs libs together: Cinit, F77init, ''
        ins = []
        for i in libs:
            ins.append(
                [l.replace("blacs", "blacs" + n) for l, n in product([i], ["Cinit", "F77init", ""])]
            )
        libs = ins + libs
        return self.check_math_libs(libs, embed)

    def _optional():
        """Check optional dependencies"""
        self.check_math_libs(OPTIONAL_DEPS, embed, optional=True)

    # parallel
    if self.get_define("ASTER_HAVE_MPI") and opts.enable_mumps:
        self.env.stash()
        try:
            # try first without blacs since now embedded by scalapack
            _scalapack()
            _optional()
            self.check_math_libs_call()
        except:
            try:
                _blacs() and _scalapack()
                _optional()
                self.check_math_libs_call()
            except:
                self.env.revert()
                _scalapack() and _blacs()
                _optional()
                self.check_math_libs_call()

    self.start_msg("Detected math libraries")
    self.end_msg(self.env[varlib])
    if self.get_define("ASTER_HAVE_MPI") and embed:
        msg = (
            "WARNING:\n"
            "    Static link with MPI libraries is not recommended.\n"
            "    Remove the option --embed-maths in case of link error.\n"
            "    See http://www.open-mpi.org/faq/?category=mpi-apps#static-mpi-apps"
        )
        Logs.warn(msg)
    if "openblas" in self.env[varlib]:
        self.define("ASTER_HAVE_OPENBLAS", 1)


@Configure.conf
def check_math_libs(self, libs, embed, optional=False):
    """Check for the first library available from 'libs'."""
    check_maths = partial(self.check_cc, uselib_store="MATH", use="MATH MPI", mandatory=False)
    if embed:

        def check_lib(lib):
            return check_maths(stlib=lib)

    else:

        def check_lib(lib):
            return check_maths(lib=lib)

    found = None
    for lib in libs:
        self.start_msg("Checking for library %s" % lib)
        if check_lib(lib=lib):
            self.end_msg("yes")
            found = lib
            break
        self.end_msg("no", color="YELLOW")
    else:
        if not optional:
            self.fatal("None of these libraries were found: %s" % libs)
    return found


@Configure.conf
def check_number_cores(self):
    """Check for the number of available cores."""
    self.start_msg("Checking for number of cores")
    try:
        self.find_program("nproc")
        try:
            res = self.cmd_and_log(["nproc"])
            nproc = int(res)
        except Errors.WafError:
            raise Errors.ConfigurationError
    except Errors.ConfigurationError:
        nproc = 1
    self.end_msg(nproc)
    self.env["NPROC"] = nproc


@Configure.conf
def get_mathlib_from_numpy(self):
    """The idea is that numpy shall be linked to blas and lapack.
    So we will try to get then using ldd if available"""
    libblas = []
    liblapack = []
    if self.env.ASTER_PLATFORM_MINGW:
        # to be rewrite for windows if needed
        return libblas, liblapack

    # numpy already checked
    cmd = self.env.PYTHON + [
        "-c",
        "\nfrom numpy.linalg import lapack_lite\nprint(lapack_lite.__file__)",
    ]
    pymodule_path = self.cmd_and_log(cmd).strip()

    self.find_program("ldd")
    ldd_env = {"LD_LIBRARY_PATH": ":".join(self.env.LIBPATH)}
    cmd = self.env.LDD + [pymodule_path]
    out = Popen(cmd, stdout=PIPE, env=ldd_env).communicate()[0].decode()

    for line in out.split("\n"):
        lib = _detect_libnames_in_ldd_line(line, LAPACK)
        if lib:
            liblapack.append(lib)
            continue
        lib = _detect_libnames_in_ldd_line(line, BLAS)
        if lib:
            libblas.append(lib)
    return libblas, liblapack


def _detect_libnames_in_ldd_line(line, libnames):
    if not list(filter(line.__contains__, libnames)):
        return None
    lib = line.split()[0].split(".", 1)[0]
    return lib[3:]


@Configure.conf
def check_math_libs_call(self, color="RED"):
    """Compile and check programs with blas/lapack, blacs and openmp"""
    self.check_math_libs_call_blas_lapack(color)
    self.check_math_libs_call_blacs(color)
    self.check_math_libs_call_openmp(color)


@Configure.conf
def check_math_libs_call_blas_lapack(self, color="RED"):
    """Compile and run a small blas/lapack program"""
    self.start_msg("Checking for a program using blas/lapack")
    try:
        ret = self.check_fc(
            fragment=blas_lapack_fragment,
            use="MPI OPENMP MATH",
            mandatory=False,
            execute=True,
            define_ret=True,
        )
        values = map(float, ret and ret.split() or [])
        ref = [10.0, 5.0]
        if list(values) != ref:
            raise Errors.ConfigurationError("invalid result: %r (expecting %r)" % (values, ref))
    except Exception as exc:
        # the message must be closed
        self.end_msg("no", color=color)
        raise Errors.ConfigurationError(str(exc))
    else:
        self.end_msg("yes")


@Configure.conf
def check_math_libs_call_blacs(self, color="RED"):
    """Compile and run a minimal blacs program"""
    if self.get_define("ASTER_HAVE_MPI"):
        self.start_msg("Checking for a program using blacs")
        try:
            self.check_fc(fragment=blacs_fragment, use="MPI OPENMP MATH", mandatory=True)
        except Exception as exc:
            # the message must be closed
            self.end_msg("no", color=color)
            raise Errors.ConfigurationError(str(exc))
        else:
            self.end_msg("yes")


@Configure.conf
def check_math_libs_call_openmp(self, color="RED"):
    """Compile and run a minimal openmp program"""
    self.start_msg("Checking for a program using omp thread")
    try:
        ret = self.check_fc(
            fragment=omp_thread_fragment,
            use="MATH OPENMP MPI",
            mandatory=True,
            execute=True,
            define_ret=True,
        )
        nbThreads = int((ret and ret.split() or [])[-1])
        refe = min(self.env["NPROC"], 2) if self.env.BUILD_OPENMP else 1
        if nbThreads < refe:
            raise ValueError("expected at least {0} thread(s)".format(nbThreads))
    except Exception as exc:
        # the message must be closed
        self.end_msg("no", color=color)
        raise Errors.ConfigurationError(str(exc))
    else:
        self.end_msg("yes (on {0} threads)".format(nbThreads))


# program testing a blas and a lapack call, output is 10.0 and 5.0
blas_lapack_fragment = r"""
subroutine test(res, res2)
    implicit none
    real(kind=8), intent(out) :: res, res2
!
    real(kind=8) :: ddot, dlapy2
    real(kind=8) :: a1(2), a2(2)
    integer  i
!
    do i = 1, 2
        a1(i) = 1.d0 * i
        a2(i) = 2.d0 * i
    end do
    res = ddot(2, a1, 1, a2,1)
    res2 = dlapy2(3.d0, 4.d0)
end subroutine test

program main
    real(kind=8) :: a, b
    call test(a, b)
    print *, a
    print *, b
end program main
"""

# program testing a blacs call, output is 0 and 1
blacs_fragment = r"""
program test_blacs
    integer iam, nprocs
    call blacs_pinfo (iam, nprocs)
    print *,iam
    print *,nprocs
end program test_blacs
"""

# program testing openmp theads
omp_thread_fragment = r"""
program hello
!$ use omp_lib
    integer total, thid
    total = 1
    thid = 0
!$omp parallel private(thid) shared(total)
!$ total = omp_get_num_threads()
!$ thid = omp_get_thread_num()
    print *, "Thread ", thid, "of ", total, "childs"
!$omp end parallel
    print *, total
end program hello
"""
