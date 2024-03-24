# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
import pathlib
import re
from dataclasses import dataclass

from waflib import Configure, Logs, TaskGen, Utils
from waflib.Context import Context
from waflib.Task import CRASHED, Task
from waflib.Tools import c, ccroot, cxx, fc


def sig_explicit_deps(self):
    """Hash `inputs` and `dep_nodes` signatures."""
    lst = []
    for node in self.inputs + self.dep_nodes:
        st = os.stat(node.abspath())
        lst.append(st.st_mtime)
        lst.append(st.st_size)
    self.m.update(Utils.h_list(lst))


# overload method for all tasks
Task.sig_explicit_deps = sig_explicit_deps


def _inputs_changed(self):
    """Tell if inputs changed."""
    for x in self.inputs + self.dep_nodes:
        for y in self.outputs:
            try:
                if os.stat(x.abspath()).st_mtime > os.stat(y.abspath()).st_mtime:
                    return True
            except:
                return True
    return False


fc_signature_native = fc.fc.signature


def signature(self):
    """By-pass signature computation if inputs haven't changed."""
    try:
        return self.cache_sig
    except AttributeError:
        pass

    if getattr(fc.fc, "_use_custom_sig", None) and not _inputs_changed(self):
        # do not compute sig_implicit_deps (and avoids scan)
        self.cache_sig = self.generator.bld.task_sigs[self.uid()]
        return self.cache_sig

    return fc_signature_native(self)


fc.fc.signature = signature

###############################################################################
# original run_str command line is store as hcode
for feature in ("program", "shlib"):
    ccroot.USELIB_VARS["fc" + feature].add("FCLINKFLAGS")
    ccroot.USELIB_VARS["c" + feature].add("CCLINKFLAGS")
    ccroot.USELIB_VARS["cxx" + feature].add("CXXLINKFLAGS")


class fcprogram(fc.fcprogram):
    """Modifying cxxprogram. Add FCLINKFLAGS."""

    run_str = "${LINK_FC} ${FCLINKFLAGS} ${LINKFLAGS} ${FCLNK_SRC_F}${SRC} ${FCLNK_TGT_F}${TGT[0].abspath()} ${RPATH_ST:RPATH} ${FCSTLIB_MARKER} ${FCSTLIBPATH_ST:STLIBPATH} ${FCSTLIB_ST:STLIB} ${FCSHLIB_MARKER} ${FCLIBPATH_ST:LIBPATH} ${FCLIB_ST:LIB} ${LDFLAGS}"


class fcshlib(fcprogram):
    """Modifying fcshlib"""

    inst_to = "${LIBDIR}"


class cprogram(c.cprogram):
    """Modifying cxxprogram. Add CCLINKFLAGS."""

    run_str = "${LINK_CC} ${CCLINKFLAGS} ${LINKFLAGS} ${CCLNK_SRC_F}${SRC} ${CCLNK_TGT_F}${TGT[0].abspath()} ${RPATH_ST:RPATH} ${FRAMEWORKPATH_ST:FRAMEWORKPATH} ${FRAMEWORK_ST:FRAMEWORK} ${ARCH_ST:ARCH} ${STLIB_MARKER} ${STLIBPATH_ST:STLIBPATH} ${STLIB_ST:STLIB} ${SHLIB_MARKER} ${LIBPATH_ST:LIBPATH} ${LIB_ST:LIB} ${LDFLAGS}"


class cshlib(cprogram):
    """Modifying cshlib"""

    inst_to = "${LIBDIR}"


class cxxprogram(cxx.cxxprogram):
    """Modifying cxxprogram. Add CXXLINKFLAGS."""

    run_str = "${LINK_CXX} ${CXXLINKFLAGS} ${LINKFLAGS} ${CXXLNK_SRC_F}${SRC} ${CXXLNK_TGT_F}${TGT[0].abspath()} ${RPATH_ST:RPATH} ${FRAMEWORKPATH_ST:FRAMEWORKPATH} ${FRAMEWORK_ST:FRAMEWORK} ${ARCH_ST:ARCH} ${STLIB_MARKER} ${STLIBPATH_ST:STLIBPATH} ${STLIB_ST:STLIB} ${SHLIB_MARKER} ${LIBPATH_ST:LIBPATH} ${LIB_ST:LIB} ${LDFLAGS}"


class cxxshlib(cxxprogram):
    """Modifying cxxshlib"""

    inst_to = "${LIBDIR}"


###############################################################################
def customize_configure_output():
    """Customize the output of configure"""

    def start_msg40(self, *k, **kw):
        """Force output on 40 columns. See :py:meth:`waflib.Context.Context.msg`"""
        if kw.get("quiet", None):
            return

        msg = kw.get("msg", None) or k[0]
        try:
            if self.in_msg:
                self.in_msg += 1
                return
        except AttributeError:
            self.in_msg = 0
        self.in_msg += 1

        self.line_just = 40  # <--- here is the change
        for x in (self.line_just * "-", msg):
            self.to_log(x)
        Logs.pprint("NORMAL", "%s :" % msg.ljust(self.line_just), sep="")

    Context.start_msg = start_msg40


customize_configure_output()


def format_error(self):
    """Write task details into a file. Print only the first line in console.
    See :py:meth:`waflib.Task.Task.format_error`"""
    text = Task.format_error(self)
    if self.hasrun == CRASHED:
        msg = getattr(self, "last_cmd", "")
        name = getattr(self.generator, "name", "")
        bldlog = osp.join(self.generator.bld.path.get_bld().abspath(), "%s.log" % name)
        try:
            os.makedirs(osp.dirname(bldlog))
        except:
            pass
        slog = ""
        try:
            open(bldlog, "w").write("task: %r\nlast command:\n%r\n" % (self, msg))
        except (OSError, IOError) as exc:
            slog = "\ncan not write the log file: %s" % str(exc)
        text = text.splitlines()[0] + "\n    task details in: {0}{1}".format(bldlog, slog)
    return text


fcprogram.format_error = format_error
cprogram.format_error = format_error
cxxprogram.format_error = format_error

# limit install lines
fun_orig = Logs.info


class CustomInfo:
    """Wrapper around `Logs.info` to minimize the number of printed lines
    during installation of Python, header and test files.
    Prints only one line for (directory, extension).
    """

    def __init__(self, prefix):
        self._dirs = set()
        self._lock = Utils.threading.Lock()
        self.prefix = prefix
        self.grouped = (
            "share/aster/tests_data",
            "share/aster/tests",
            "share/locale/aster",
            "lib/aster/code_aster",
            "lib/aster/run_aster",
        )

    def __call__(self, *args, **kwargs):
        group = False
        if len(args) == 6 and ("+ install" in args[0] or "- install" in args[0]):
            group = True
        if not group:
            fun_orig(*args, **kwargs)
            return

        dst = args[3]
        dirn = osp.dirname(dst)
        ext = osp.splitext(dst)[-1]
        key = dirn, ext
        grouped, grp = self._grouped(dirn)
        if grouped:
            key = grp
            ext = ""
        self._lock.acquire()
        show = key not in self._dirs
        self._dirs.add(key)
        self._lock.release()
        if show:
            args = list(args)
            args[3] = osp.join(grp, "*" + ext)
            args[5] = osp.dirname(args[5])
            args[0] = args[0].replace(" (from %s)", "")
            args.pop()
            fun_orig(*args, **kwargs)

    def _grouped(self, dirn):
        for i in self.grouped:
            if i in dirn:
                return True, osp.join(self.prefix, i)
        return False, dirn


def build(self):
    if Logs.verbose < 1:
        Logs.info = CustomInfo(self.env.PREFIX)


# support for the "dynamic_source" attribute
@TaskGen.feature("c", "cxx")
@TaskGen.before("process_source", "process_rule")
def dynamic_post(self):
    """
    bld(dynamic_source='*.c', ...)
        will search for source files to add to the attribute 'source'.

    bld(dynamic_source='*.c', dynamic_incpaths='include', ...)
        will search for 'include' in the parent of every new source and
        add it in INCLUDES paths.
    """
    if not getattr(self, "dynamic_source", None):
        return
    self.source = Utils.to_list(self.source)
    get_srcs = self.path.get_bld().ant_glob
    added = get_srcs(self.dynamic_source, remove=False, quiet=True)
    self.source.extend(added)
    for node in added:
        node.sig = Utils.h_file(node.abspath())
        if getattr(self, "dynamic_incpaths", None):
            incpath = node.parent.find_node(self.dynamic_incpaths)
            if incpath:
                incpath.sig = incpath.abspath()
                self.env.append_value("INCLUDES", [incpath.abspath()])
                incs = incpath.get_bld().ant_glob("**/*.h*", quiet=True)
                for node in incs:
                    node.sig = Utils.h_file(node.abspath())


###############################################################################
@Configure.conf
def safe_remove(self, var, value):
    """Remove 'value' from the variable, remove duplicates"""
    if type(self.env[var]) is not list:
        return
    self.env[var] = self.remove_duplicates(self.env[var])
    while value in self.env[var]:
        self.env[var].remove(value)


@Configure.conf
def remove_flags(self, var, flags):
    """Remove `flags` from `env[var]`."""
    if not isinstance(self.env[var], list):
        return
    regexps = [re.compile(re.escape(flag) + ".*") for flag in flags]
    # TODO itertools.chain?
    flatten = lambda l: [item for sublist in l for item in sublist]
    fl = []
    for regexp in regexps:
        fl.extend(flatten(filter(None, map(regexp.findall, self.env[var]))))
    self.env[var] = list(set([i for i in self.env[var] if i not in fl]))


@Configure.conf
def remove_optflags(self, type_flags):
    """Remove optimisation flags from the `type_flags`/* variables,
    remove duplicates"""
    for var in self.env:
        if var.startswith(type_flags):
            if not isinstance(self.env[var], (list, tuple)):
                self.env[var] = [self.env[var]]
            self.env[var] = self.remove_duplicates(self.env[var])
            self.env[var] = [i for i in self.env[var] if not i.startswith("-O")]


@Configure.conf
def remove_duplicates(self, list_in):
    """Return the list by removing the duplicated elements
    and by keeping the order. It ignores empty values."""
    dset = set()
    # relies on the fact that dset.add() always returns None.
    return [path for path in list_in if path not in dset and not dset.add(path)]


@TaskGen.feature("cshlib")
@TaskGen.before_method('apply_link')
def fix_bibc_deps(self):
    for task in self.bld.get_tgen_by_name("asterbibc").tasks:
        Logs.info(f"{task.__class__.__name__=}")
        Logs.info(f"{task.inputs=}")
        Logs.info(f"{task=}")


@dataclass
class ExtraDep:
    source_task: Task
    dep_nodes: list[Task]


def fix_specific_c_bibfor_includes(self):
    added_links = {
        "cshlib": ["bibfor/utilitai/utgtme.F90", "bibfor/utilitai/utptme.F90", "bibfor/parallel/asmpi_warn.F90",
                   "bibfor/utilitai/fclose.F90", "bibfor/utilitai/ulopen.F90", "bibfor/supervis/affich.F90",
                   "bibfor/supervis/impers.F90", "bibfor/prepost/mdnoma.F90", "bibfor/prepost/mdnoch.F90",
                   "bibfor/modelisa/rcvale_wrap.F90", "bibfor/modelisa/gmardm.F90", "bibfor/utilifor/postkutil.F90",
                   "bibfor/jeveux/jemarq.F90", "bibfor/supervis/getcon.F90", "bibfor/jeveux/jedema.F90",
                   "bibfor/jeveux/jedetr.F90", "bibfor/supervis/putcon.F90", "bibfor/supervis/tailsd.F90",
                   "bibfor/resu_util/rsacch.F90", "bibfor/resu_util/rsacpa.F90", "bibfor/jeveux/jelst3.F90",
                   "bibfor/jeveux/jelira.F90", "bibfor/jeveux/jeexin.F90", "bibfor/supervis/gcncon.F90",
                   "bibfor/utilifor/abortf.F90", "bibfor/supervis/sigcpu.F90", "bibfor/utilitai/utmfpe.F90",
                   "bibfor/utilifor/utmess_cwrap.F90", "bibfor/utilifor/utmess_core.F90"],

        "cxxshlib": ["bibfor/jeveux/wkvectc.F90", "bibfor/jeveux/jeecra_wrap.F90", "bibfor/jeveux/juveca.F90",
                     "bibc/supervis/aster_utils.c", "bibfor/jeveux/jedema.F90", "bibfor/jeveux/jemarq.F90",
                     "bibfor/jeveux/jexnom.F90", "bibfor/matrix/delete_matrix.F90"]
    }
    top_dir = pathlib.Path(self.bld.top_dir)
    lib_deps = {}
    all_outputs = []
    for task in self.bld.get_tgen_by_name("asterbibfor").tasks:
        if task.__class__.__name__ != 'fc':
            continue
        fp_in = pathlib.Path(task.inputs[0].abspath())
        lib_deps[fp_in] = task
        all_outputs.extend(task.outputs)

    for task in self.bld.get_tgen_by_name("asterbibc").tasks:
        tname = task.__class__.__name__
        if tname != 'cshlib':
            continue
        # task.inputs.extend(all_outputs)
        deps = added_links[tname]
        for dep in deps:
            fp_in = pathlib.Path(task.inputs[0].abspath())
            dep_fp = top_dir / dep
            if dep_fp not in lib_deps:
                Logs.error(f"Dependency {dep_fp} not found in bibfor_deps")
                continue
            dep_task = lib_deps[dep_fp]
            lib_deps[fp_in] = task
            task.inputs.extend(dep_task.outputs)
            Logs.info(f"Adding {dep_task} to {task}")
            Logs.info(f"{task.inputs=}")

    for task in self.bld.get_tgen_by_name("asterbibcxx").tasks:
        tname = task.__class__.__name__
        if tname != 'cxxshlib':
            continue
        # task.inputs.extend(all_outputs)
        deps = added_links[tname]
        for dep in deps:
            dep_fp = top_dir / dep
            if dep_fp not in lib_deps:
                Logs.error(f"Dependency {dep_fp} not found in bibfor_deps")
                continue
            dep_task = lib_deps[dep_fp]
            task.inputs.extend(dep_task.outputs)
            Logs.info(f"Adding {dep_task} to {task}")
            Logs.info(f"{task.inputs=}")


@TaskGen.feature("cshlib")
@TaskGen.after_method('apply_link')
def def_prep_fc_c_linking(self):
    fix_specific_c_bibfor_includes(self)

    fc_task = None
    for task in self.bld.get_tgen_by_name("asterbibfor").tasks:
        if task.__class__.__name__ != 'fcshlib':
            continue
        Logs.info(f"{task.__class__.__name__=}")
        Logs.info(f"{task.outputs=}")
        # Logs.info(f"{task.inputs=}")
        # Logs.info(f"{task=}")
        fc_task = task

    c_task = None
    for task in self.bld.get_tgen_by_name("asterbibc").tasks:
        if task.__class__.__name__ != 'cshlib':
            continue
        Logs.info(f"{task.__class__.__name__=}")
        Logs.info(f"{task.outputs=}")
        # Logs.info(f"{task.inputs=}")
        # Logs.info(f"{task=}")
        c_task = task

    cxx_task = None
    for task in self.bld.get_tgen_by_name("asterbibcxx").tasks:
        if task.__class__.__name__ != 'cxxshlib':
            continue
        Logs.info(f"{task.__class__.__name__=}")
        Logs.info(f"{task.outputs=}")
        # Logs.info(f"{task.inputs=}")
        # Logs.info(f"{task=}")
        cxx_task = task

    # Print depends order for each task
    Logs.info(f"{fc_task.dep_nodes=}")
    Logs.info(f"{c_task.dep_nodes=}")
    Logs.info(f"{cxx_task.dep_nodes=}")

    # These will Force linking order: bibfor -> bibc -> bibcxx
    c_task.dep_nodes.append(fc_task)
    cxx_task.dep_nodes.append(fc_task)
    cxx_task.dep_nodes.append(c_task)

    # Not quite sure about these two lines
    # c_task.inputs.extend(fc_task.outputs)
    # cxx_task.inputs.extend(fc_task.outputs)

    Logs.info(f"{c_task.priority()=}")
    Logs.info(f"{cxx_task.priority()=}")
    Logs.info(f"{fc_task.priority()=}")
