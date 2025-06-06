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
import pathlib
import platform
from pathlib import PureWindowsPath
from subprocess import PIPE, Popen

from waflib import Configure, Errors, Logs


def options(self):
    self.load("python")  # --nopyc/--nopyo are enabled below


def configure(self):
    self.check_python()
    self.check_numpy()
    self.check_asrun()
    self.check_mpi4py()


###############################################################################
@Configure.conf
def check_python(self):
    self.load("python")
    self.check_python_version((3, 5, 0))
    if platform.system() == "Windows":
        path = self.env["PATH"]
        include_dir = pathlib.Path(self.env.PREFIX) / "include"
        self.env["PATH"] = f"{path};{include_dir.as_posix()}"
    else:
        self.check_python_headers()

    if "icc" in self.env.CC_NAME.lower():
        self.env["LIB_PYEXT"] = list(set(self.env["LIB_PYEXT"]))
        # Best is to clear PYEMBED and PYEXT {c/cxx}flags
        for lang in ("CFLAGS", "CXXFLAGS"):
            for feat in ("PYEMBED", "PYEXT"):
                self.env[lang + "_" + feat] = []
    cfgext = self.env["CONFIG_PARAMETERS"].get("cfgext", "yaml")
    if cfgext in ("", "yaml"):
        try:
            self.check_python_module("yaml")
        except Errors.ConfigurationError:
            cfgext = "json"
    self.env["CFG_EXT"] = cfgext


@Configure.conf
def check_numpy(self):
    if not self.env["PYTHON"]:
        self.fatal("load python tool first")
    self.check_python_module("numpy")
    self.check_numpy_headers()


@Configure.conf
def check_numpy_headers(self):
    if not self.env["PYTHON"]:
        self.fatal("load python tool first")
    self.start_msg("Checking for numpy includes")
    # retrieve includes dir from numpy module
    cmd = self.env.PYTHON + ["-c", "\nimport numpy\nprint(numpy.get_include())"]
    numpy_includes = self.cmd_and_log(cmd, shell=False).strip()
    self.end_msg(numpy_includes)
    self.start_msg("Checking for numpy arrayobject.h")

    extra_flags = dict()
    if platform.system() == "Windows":
        library_prefix_ = pathlib.Path(self.env.PREFIX)
        prefix_ = library_prefix_.parent
        python_include_dir = prefix_ / "include"
        python_libs_dir = prefix_ / "libs"
        numpy_includes.append(python_include_dir.as_posix())
        numpy_includes.append(python_libs_dir.as_posix())
        extra_flags.update(dict(linkflags=[f"/LIBPATH:{python_libs_dir}", f"/LIBPATH:{python_include_dir}"]))
        Logs.info(f"MSVC {numpy_includes=}")

    if self.is_defined("ASTER_PLATFORM_MINGW"):
        incs = PureWindowsPath(numpy_includes)
        parts = list(incs.parts)
        if incs.anchor:
            parts[0] = incs.root
        for i, sub in enumerate(parts):
            if sub == "lib":
                parts[i] = "Lib"
        numpy_includes = PureWindowsPath(*parts).as_posix()
    # check the given includes dirs
    self.check(
        feature="c",
        header_name="Python.h numpy/arrayobject.h",
        includes=numpy_includes,
        use=["PYEXT"],
        uselib_store="NUMPY",
        errmsg="Could not find the numpy development headers",
        **extra_flags
    )
    self.end_msg(numpy_includes)


@Configure.conf
def check_asrun(self):
    if not self.env["PYTHON"]:
        self.fatal("load python tool first")
    try:
        self.check_python_module("asrun")
    except Errors.WafError:
        # optional
        pass


@Configure.conf
def check_mpi4py(self):
    if not self.env.BUILD_MPI:
        return
    if not self.env["PYTHON"]:
        self.fatal("load python tool first")
    try:
        self.check_python_module("mpi4py")
    except Errors.ConfigurationError:
        if self.options.with_py_mpi4py:
            raise


@Configure.conf
def check_optimization_python(self):
    self.setenv("debug")
    self.env["PYC"] = self.env["PYO"] = 0
    self.setenv("release")
    self.env["PYC"] = self.env["PYO"] = 0


def _get_default_pythonpath(python):
    """Default sys.path should be added into PYTHONPATH"""
    env = os.environ.copy()
    env["PYTHONPATH"] = ""
    proc = Popen([python, "-c", "import sys; print(sys.path)"], stdout=PIPE, env=env)
    system_path = eval(proc.communicate()[0])
    return system_path
