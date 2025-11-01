import os
import pathlib
import platform
import shutil
from dataclasses import dataclass

from waflib import Logs, TaskGen, Task, Errors

from waftools.clangdb import build_clang_compilation_db


class msvclibgen(Task.Task):
    # Use a minimal run_str; we will compose the full command with a response file in exec_command
    run_str = "LIB.exe /OUT:${TGT}"
    color = "BLUE"
    before = ["cshlib", "cxxshlib", "fcshlib"]
    use_msvc_entry = False

    def __str__(self):
        # Provide a concise task banner to avoid dumping thousands of inputs
        try:
            out_name = self.outputs[0].name if self.outputs else "<no-output>"
            n_inputs = len(self.inputs) if hasattr(self, "inputs") else 0
            return f"msvclibgen: {out_name} ({n_inputs} objects)"
        except Exception:
            # Fall back to default if anything unexpected happens
            return super().__str__()

    def exec_command(self, cmd, **kw):
        """Execute the command with a response file to avoid massive command lines/logs."""
        output_fp = pathlib.Path(self.outputs[0].abspath())
        output_fp.parent.mkdir(parents=True, exist_ok=True)
        obld = self.generator.bld
        root_path = pathlib.Path(obld.root.abspath()).resolve().absolute()
        clean_name_map = {
            "bibc": self.env.BIBC_DEF,
            "bibcxx": self.env.BIBCXX_DEF,
            "bibfor": self.env.BIBFOR_DEF,
            "bibfor_ext": self.env.BIBFOR_EXT_DEF,
            "AsterGC": self.env.ASTERGC_DEF,
            "aster": self.env.ASTER_DEF,
            "mfront": self.env.MFRONT_DEF,
        }
        clean_name = output_fp.stem.replace("_gen", "")
        # libs directory for dependent libs
        libs_dir = pathlib.Path(self.env.PREFIX).resolve().absolute().parent / "libs"
        # Base options
        opts = ["/NOLOGO", "/MACHINE:X64", "/SUBSYSTEM:CONSOLE", f"/LIBPATH:{libs_dir}"]
        # DEF file
        def_file = root_path / clean_name / f"{clean_name}.def"
        if clean_name.endswith("proxy"):
            def_file = root_path / "msvc/c_entrypoints" / f"{clean_name}.def"
            opts += [f"/DEF:{def_file}"]
        else:
            def_file_v = clean_name_map.get(clean_name, None)
            if def_file_v is not None:
                def_file = root_path / def_file_v
                Logs.debug(f"Using def file {def_file=}")
            opts += [f"/DEF:{def_file}"]

        # Build a response file containing all input object files
        # self.inputs are Nodes (previous compile tasks' outputs)
        src_paths = [str(n.abspath()) for n in self.inputs]
        rsp_path = output_fp.with_suffix(".rsp")
        with open(rsp_path, "w", encoding="utf-8") as rsp:
            # one entry per line keeps it readable
            for p in src_paths:
                rsp.write(p)
                rsp.write("\n")

        # Compose compact command: LIB.exe /OUT:<out> <opts> @<rsp>
        cmd_compact = [cmd[0], cmd[1], *opts, f"@{rsp_path}"]

        # Log a concise summary at info level
        Logs.info(
            f"Generating {output_fp.name} with {len(src_paths)} objects (using response file '{rsp_path}')."
        )
        Logs.debug(f"Response file: {rsp_path}")

        # Execute
        ret = super().exec_command(cmd_compact, **kw)
        return ret


@dataclass
class TaskObject:
    libtask: Task
    program_task: Task
    tasks: list[Task]
    task_gen: TaskGen


@dataclass
class LibTask:
    asterbibc: TaskObject
    asterbibcxx: TaskObject
    asterbibfor: TaskObject
    asterbibfor_ext: TaskObject
    astergc: TaskObject
    asterlib: TaskObject

    def all_tasks_ready(self) -> bool:
        for key, task in self.__dict__.items():
            if len(task.tasks) == 0:
                Logs.error(f"No tasks found for {task.task_gen.get_name()=}")
                return False
        return True

    def get_missing_tasks(self) -> list[str]:
        missing_tasks = []
        for key, task in self.__dict__.items():
            if len(task.tasks) == 0:
                missing_tasks.append(task.task_gen.get_name())

        return missing_tasks


def get_task_object(bld: TaskGen, taskgen_name: str, compiler_prefix: str) -> TaskObject:
    lib_task = None
    program_task = None
    tasks = []
    task_gen = bld.get_tgen_by_name(taskgen_name)
    for task in task_gen.tasks:
        if task.__class__.__name__ == f"{compiler_prefix}shlib":
            lib_task = task
        if task.__class__.__name__ == f"{compiler_prefix}program":
            program_task = task
        if task.__class__.__name__ == compiler_prefix:
            tasks.append(task)

    if len(tasks) == 0:
        Errors.WafError(f"No tasks found for {taskgen_name}")

    return TaskObject(lib_task, program_task, tasks, task_gen)


def extract_main_tasks(self: TaskGen.task_gen) -> LibTask:
    Logs.debug(f"Extracting main tasks for {self.get_name()=}")

    bld = self.bld
    c_task_object = get_task_object(bld, "asterbibc", "c")
    cxx_task_object = get_task_object(bld, "asterbibcxx", "cxx")
    fc_task_object = get_task_object(bld, "asterbibfor", "fc")
    fc_ext_task_object = get_task_object(bld, "asterbibfor_ext", "fc")
    gc_task_object = get_task_object(bld, "astergc", "cxx")
    c_aster_object = get_task_object(bld, "asterlib", "cxx")

    return LibTask(c_task_object, cxx_task_object, fc_task_object, fc_ext_task_object, gc_task_object, c_aster_object)


def create_msvclibgen_task(self, lib_name: str, input_tasks) -> Task:
    # Create a task for MSVC lib generation for C
    # the task takes in the outputs of all C tasks
    # it's outputs are bibc.lib and bibc.exp located in the build directory
    bld_path = pathlib.Path(self.bld.bldnode.abspath()).resolve().absolute()
    if lib_name == "aster":
        lib_output_file_path = bld_path / "bibc" / "aster.lib"
    elif lib_name == "bibfor_ext":
        lib_output_file_path = bld_path / "bibfor" / "bibfor_ext.lib"
    elif lib_name == "AsterGC":
        lib_output_file_path = bld_path / "libs" / "AsterGC.lib"
    elif lib_name.endswith("proxy"):
        Logs.debug(f"input_tasks: {input_tasks=}")
        lib_output_file_path = bld_path / "msvc" / f"{lib_name}.lib"
    else:
        lib_output_file_path = bld_path / lib_name / f"{lib_name}.lib"

    # create nodes for the output files
    bib_lib_output_file_node = self.bld.bldnode.make_node(lib_output_file_path.relative_to(bld_path).as_posix())

    msvc_libgen_task = self.create_task("msvclibgen")
    msvc_libgen_task.inputs = input_tasks
    msvc_libgen_task.env = self.env
    msvc_libgen_task.dep_nodes = input_tasks
    msvc_libgen_task.outputs = [bib_lib_output_file_node]

    Logs.debug(f"{msvc_libgen_task.outputs=}")

    return msvc_libgen_task


def run_mvsc_lib_gen(self, task_obj: LibTask):
    Logs.info("Generating MSVC import libraries")
    clib_task = task_obj.asterbibc.libtask
    cxxlib_task = task_obj.asterbibcxx.libtask
    fclib_task = task_obj.asterbibfor.libtask
    fcext_lib_task = task_obj.asterbibfor_ext.libtask
    gc_lib_task = task_obj.astergc.libtask
    aster_task = task_obj.asterlib.libtask

    Logs.debug(f"Before removal: {clib_task.outputs=}")

    # Lib files are created by MSVC lib generation, so will remove these from the shlib outputs
    clib_task.outputs = [o for o in clib_task.outputs if o.suffix() != ".lib"]
    cxxlib_task.outputs = [o for o in cxxlib_task.outputs if o.suffix() != ".lib"]
    fclib_task.outputs = [o for o in fclib_task.outputs if o.suffix() != ".lib"]
    fcext_lib_task.outputs = [o for o in fcext_lib_task.outputs if o.suffix() != ".lib"]
    gc_lib_task.outputs = [o for o in gc_lib_task.outputs if o.suffix() != ".lib"]
    aster_task.outputs = [o for o in aster_task.outputs if o.suffix() != ".lib"]

    Logs.debug(f"After removal: {clib_task.outputs=}")
    Logs.debug(f"After removal: {fclib_task.outputs=}")
    Logs.debug(f"After removal: {fcext_lib_task.outputs=}")
    Logs.debug(f"After removal: {gc_lib_task.outputs=}")
    Logs.debug(f"After removal: {cxxlib_task.outputs=}")
    Logs.debug(f"After removal: {aster_task.outputs=}")

    c_input_tasks = [ctask.outputs[0] for ctask in task_obj.asterbibc.tasks]
    cxx_input_tasks = [cxxtask.outputs[0] for cxxtask in task_obj.asterbibcxx.tasks]
    fc_input_tasks = [fctask.outputs[0] for fctask in task_obj.asterbibfor.tasks]
    fc_ext_input_tasks = [fctask.outputs[0] for fctask in task_obj.asterbibfor_ext.tasks]
    gc_input_tasks = [gctask.outputs[0] for gctask in task_obj.astergc.tasks]
    aster_input_tasks = [ctask.outputs[0] for ctask in task_obj.asterlib.tasks if ctask.outputs[0].suffix() == ".o"]
    Logs.debug(f"{aster_input_tasks=}")

    if len(aster_input_tasks) == 0:
        raise Errors.WafError("Failed MSVC lib generation: No aster input tasks found")

    clib_lib_task = create_msvclibgen_task(self, "bibc", c_input_tasks)
    bibcxx_lib_task = create_msvclibgen_task(self, "bibcxx", cxx_input_tasks)
    bibfor_lib_task = create_msvclibgen_task(self, "bibfor", fc_input_tasks)
    bibfor_ext_lib_task = create_msvclibgen_task(self, "bibfor_ext", fc_ext_input_tasks)
    gc_lib_task_gen = create_msvclibgen_task(self, "AsterGC", gc_input_tasks)
    bibaster_lib_task = create_msvclibgen_task(self, "aster", aster_input_tasks)

    Logs.debug(f"{clib_lib_task.outputs=}")
    Logs.debug(f"{fclib_task.outputs=}")
    Logs.debug(f"{fcext_lib_task.outputs=}")
    Logs.debug(f"{gc_lib_task.outputs=}")
    Logs.debug(f"{bibaster_lib_task.outputs=}")

    # filter out all non-lib files
    clib_task_outputs = [x for x in clib_lib_task.outputs if x.suffix() == ".lib"]
    fclib_task_outputs = [x for x in bibfor_lib_task.outputs if x.suffix() == ".lib"]
    fcext_lib_task_outputs = [x for x in bibfor_ext_lib_task.outputs if x.suffix() == ".lib"]
    gc_task_outputs = [x for x in gc_lib_task_gen.outputs if x.suffix() == ".lib"]
    bibaster_task_outputs = [x for x in bibaster_lib_task.outputs if x.suffix() == ".lib"]
    bibcxx_task_outputs = [x for x in bibcxx_lib_task.outputs if x.suffix() == ".lib"]

    Logs.debug(f"{clib_task_outputs=}")
    Logs.debug(f"{fclib_task_outputs=}")
    Logs.debug(f"{fcext_lib_task_outputs=}")
    Logs.debug(f"{gc_task_outputs=}")
    Logs.debug(f"{bibaster_task_outputs=}")
    Logs.debug(f"{bibcxx_task_outputs=}")

    # Set up library dependencies according to the dependency chain
    # Strategy: Generate all .lib files first (they only need object files),
    # then link all DLLs (which need the .lib files)

    # Note: msvclibgen tasks (bibfor_lib_task, etc.) only need object files as inputs.
    # They do NOT depend on other .lib files, so NO set_run_after between them.
    # This avoids circular dependencies.
    # Only the DLL link tasks need to wait for .lib generation to complete.

    # All DLL link tasks must wait for ALL .lib generation tasks to complete
    all_lib_gen_tasks = [bibfor_lib_task, bibfor_ext_lib_task, clib_lib_task,
                         bibcxx_lib_task, gc_lib_task_gen, bibaster_lib_task]

    # bibfor.dll depends on bibcxx.lib and bibc.lib
    fclib_task.inputs += bibcxx_task_outputs + clib_task_outputs
    for lib_task in all_lib_gen_tasks:
        fclib_task.set_run_after(lib_task)

    # bibfor_ext.dll depends on bibfor.lib, bibcxx.lib, and bibc.lib
    fcext_lib_task.inputs += fclib_task_outputs + bibcxx_task_outputs + clib_task_outputs
    for lib_task in all_lib_gen_tasks:
        fcext_lib_task.set_run_after(lib_task)

    # bibc.dll depends on bibfor.lib, bibfor_ext.lib, and bibcxx.lib
    clib_task.inputs += fclib_task_outputs + fcext_lib_task_outputs + bibcxx_task_outputs
    for lib_task in all_lib_gen_tasks:
        clib_task.set_run_after(lib_task)

    # AsterGC.dll depends on bibfor.lib and bibc.lib
    gc_lib_task.inputs += fclib_task_outputs + clib_task_outputs
    for lib_task in all_lib_gen_tasks:
        gc_lib_task.set_run_after(lib_task)

    # bibcxx.dll depends on aster.lib, bibc.lib, bibfor.lib, bibfor_ext.lib, and AsterGC.lib
    cxxlib_task.inputs += bibaster_task_outputs + clib_task_outputs + fclib_task_outputs + fcext_lib_task_outputs + gc_task_outputs
    for lib_task in all_lib_gen_tasks:
        cxxlib_task.set_run_after(lib_task)

    # aster.dll depends on bibfor.lib, bibfor_ext.lib, bibc.lib, bibcxx.lib, and AsterGC.lib
    aster_task.inputs += fclib_task_outputs + fcext_lib_task_outputs + clib_task_outputs + bibcxx_task_outputs + gc_task_outputs
    for lib_task in all_lib_gen_tasks:
        aster_task.set_run_after(lib_task)

    bibc_dll = [o for o in clib_task.outputs if o.suffix() == ".dll"][0]
    bibcxx_dll = [o for o in cxxlib_task.outputs if o.suffix() == ".dll"][0]
    Logs.info(f"{type(bibc_dll)=}{bibc_dll=}")
    Logs.info(f"{type(bibcxx_dll)=}{bibcxx_dll=}")

    Logs.info("Successfully ran MSVC lib generation")


_lib_task_obj: LibTask | None = None
_task_done = False
_compiler_map = {
    "asterlib": "cxx",
    "asterbibc": "c",
    "asterbibfor": "fc",
    "asterbibfor_ext": "fc",
    "astergc": "cxx",
    "asterbibcxx": "cxx"
}

_proxy_tasks = set()


@TaskGen.feature("cxxshlib")
@TaskGen.after_method("apply_link", "propagate_uselib_vars")
def make_libs_for_entrypoints(self) -> None:
    if platform.system() != "Windows":
        return
    name = self.get_name()
    if not name.endswith("proxy"):
        return
    global _proxy_tasks
    if name in _proxy_tasks:
        return
    _proxy_tasks.add(name)

    aster_object = get_task_object(self.bld, name, "cxx")
    c_input_tasks = [ctask.outputs[0] for ctask in aster_object.tasks]
    clib_task = aster_object.libtask

    Logs.debug(f"{name=},{c_input_tasks=}, {clib_task.outputs=}")

    clib_task.outputs = [o for o in clib_task.outputs if o.suffix() != ".lib"]
    Logs.debug(f"{clib_task.outputs=} after removal")

    aster_msvc_lib_task = create_msvclibgen_task(self, name, c_input_tasks)

    clib_task_outputs = [x for x in aster_msvc_lib_task.outputs if x.suffix() == ".lib"]

    clib_task.inputs += clib_task_outputs
    Logs.info(f"{clib_task.inputs=}")


@TaskGen.feature("cxxshlib", "fcshlib", "cshlib")
@TaskGen.after_method("apply_link", "propagate_uselib_vars")
def set_flags(self) -> None:
    if platform.system() != "Windows":
        return
    name = self.get_name()
    args = []
    bld_path = pathlib.Path(self.bld.bldnode.abspath()).resolve().absolute()
    if name == "asterbibc":
        archive_name = "bibc"
    elif name == "asterbibcxx":
        archive_name = "bibcxx"
    elif name == "asterbibfor":
        archive_name = "bibfor"
    elif name == "asterbibfor_ext":
        archive_name = "bibfor_ext"
    elif name == "astergc":
        archive_name = "AsterGC"
    elif name == "asterlib":
        archive_name = "aster"
    elif name.endswith("proxy"):
        archive_name = name
        conda_dir = bld_path / "msvc"
        args += [f"/LIBPATH:{conda_dir.as_posix()}"]
    else:
        Logs.info(f"Skipping {name=}")
        return None


    if name == "astergc":
        archive_dir = bld_path / 'libs'
    elif name == "asterbibfor_ext":
        archive_dir = bld_path / 'bibfor'
    else:
        archive_dir = bld_path / archive_name

    # Add LIBPATH for the library's own directory
    args += [f"/LIBPATH:{archive_dir.as_posix()}", f"/WHOLEARCHIVE:{archive_name}.lib", f"{archive_name}.exp"]

    # Add LIBPATH for all dependent libraries so the linker can find their import libraries
    # All libraries potentially depend on each other, so add all paths
    dependent_lib_paths = [
        bld_path / 'bibc',
        bld_path / 'bibcxx',
        bld_path / 'bibfor',
        bld_path / 'libs',  # for AsterGC
    ]

    for lib_path in dependent_lib_paths:
        if lib_path != archive_dir:  # Don't duplicate the library's own path
            args += [f"/LIBPATH:{lib_path.as_posix()}"]

    Logs.debug(f"{archive_name=} extra flags {args} for {name=}")
    self.link_task.env.append_unique("LINKFLAGS", args)
    # log the entire LINKFLAGS
    Logs.debug(f"{archive_name=}: {self.link_task.env.LINKFLAGS=}")

@TaskGen.feature("cxxshlib", "fcshlib", "cshlib")
@TaskGen.after_method("apply_link")
def make_msvc_modifications(self: TaskGen.task_gen):
    name = self.get_name()
    Logs.info(f"task: {name=}")

    global _lib_task_obj
    global _task_done
    if _task_done:
        Logs.info("Task already done")
        return

    if _lib_task_obj is None:
        _lib_task_obj = extract_main_tasks(self)
    lib_task_obj = _lib_task_obj

    if name in _compiler_map.keys():
        Logs.info(f"setting {name=}")
        aster_object = get_task_object(self.bld, name, _compiler_map[name])
        if name == "asterlib" and len(aster_object.tasks) == 0:
            # For some reason the asterlib compile tasks are forced to be compiled with the c compiler despite it
            # being assigned cxx compiler in the wscript
            aster_object_c = get_task_object(self.bld, name, "c")
            aster_object.tasks = aster_object_c.tasks
        setattr(lib_task_obj, name, aster_object)

    are_tasks_ready = lib_task_obj.all_tasks_ready()
    Logs.info(f"{are_tasks_ready=}")

    if are_tasks_ready is False:
        missing_task_names = lib_task_obj.get_missing_tasks()
        Logs.error(f"Missing tasks @{name} after: {missing_task_names}")
        return

    Logs.info("All tasks are ready")
    run_mvsc_lib_gen(self, lib_task_obj)
    build_clang_compilation_db(
        self,
        lib_task_obj.asterbibc.tasks,
        lib_task_obj.asterbibcxx.tasks + lib_task_obj.asterlib.tasks + lib_task_obj.astergc.tasks,
        lib_task_obj.asterbibfor.tasks + lib_task_obj.asterbibfor_ext.tasks,
    )
    _task_done = True
