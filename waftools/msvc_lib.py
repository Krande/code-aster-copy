import os
import pathlib
import platform
import shutil
from dataclasses import dataclass

from waflib import Logs, TaskGen, Task, Errors

from waftools.clangdb import build_clang_compilation_db


class msvclibgen(Task.Task):
    run_str = "LIB.exe /OUT:${TGT} ${SRC}"
    color = "BLUE"
    before = ["cshlib", "cxxshlib", "fcshlib"]

    def exec_command(self, cmd, **kw):
        """Execute the command"""
        output_fp = pathlib.Path(self.outputs[0].abspath())
        output_fp.parent.mkdir(parents=True, exist_ok=True)
        obld = self.generator.bld
        root_path = pathlib.Path(obld.root.abspath()).resolve().absolute()

        clean_name = output_fp.stem.replace("_gen", "")
        # This is a hack to copy the generated lib to the build directory
        if clean_name == "aster":
            def_file = root_path / "bibc" / "aster.def"
        else:
            def_file = root_path / clean_name / f"{clean_name}.def"

        # Location of python 3.11 libs
        libs_dir = pathlib.Path(self.env.PREFIX).resolve().absolute().parent / "libs"
        opts = ["/NOLOGO", "/MACHINE:X64", "/SUBSYSTEM:CONSOLE", f"/LIBPATH:{libs_dir}", f"/DEF:{def_file}"]

        cmd = cmd[:2] + opts + cmd[2:]
        # write a copy of the inputs to a file
        with open(output_fp.with_name(f"{output_fp.stem}_in.txt"), "w") as f:
            f.write("\n".join([str(x) for x in cmd]))

        ret = super().exec_command(cmd, **kw)
        return ret


def create_symlink(source, link_name):
    """
    Attempts to create a symbolic link on Windows.

    Args:
    source (str): The path to the target file or directory.
    link_name (str): The path where the symbolic link should be created.

    Returns:
    bool: True if the symlink was created successfully, False otherwise.
    """
    try:
        os.symlink(source, link_name)
        print(f"Symlink created successfully: {link_name} -> {source}")
        return True
    except OSError as e:
        print(f"Failed to create symlink: {e}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False


class msvc_symlink_installer(Task.Task):
    ext_in = "install"

    def run(self):
        for i, in_file in enumerate(self.inputs):
            in_file_fp = pathlib.Path(in_file.abspath())
            output_fp = pathlib.Path(self.outputs[i].abspath())
            Logs.info(f"Creating symlink: {in_file_fp} -> {output_fp}")
            result = create_symlink(in_file_fp, output_fp)
            if not result:
                shutil.copy(in_file_fp, output_fp)
                Logs.info(f"Failed to create symlink: {in_file_fp} -> {output_fp}, therefore copying file instead")
        return 0


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
    Logs.info(f"Extracting main tasks for {self.get_name()=}")

    bld = self.bld
    c_task_object = get_task_object(bld, "asterbibc", "c")
    cxx_task_object = get_task_object(bld, "asterbibcxx", "cxx")
    fc_task_object = get_task_object(bld, "asterbibfor", "fc")
    c_aster_object = get_task_object(bld, "asterlib", "cxx")

    return LibTask(c_task_object, cxx_task_object, fc_task_object, c_aster_object)


def create_msvclibgen_task(self, lib_name: str, input_tasks) -> Task:
    # Create a task for MSVC lib generation for C
    # the task takes in the outputs of all C tasks
    # it's outputs are bibc.lib and bibc.exp located in the build directory
    bld_path = pathlib.Path(self.bld.bldnode.abspath()).resolve().absolute()
    if lib_name == "aster":
        lib_output_file_path = bld_path / "bibc" / "aster.lib"
    else:
        lib_output_file_path = bld_path / lib_name / f"{lib_name}.lib"

    # create nodes for the output files
    bib_lib_output_file_node = self.bld.bldnode.make_node(lib_output_file_path.relative_to(bld_path).as_posix())

    msvc_libgen_task = self.create_task("msvclibgen")
    msvc_libgen_task.inputs = input_tasks
    msvc_libgen_task.env = self.env
    msvc_libgen_task.dep_nodes = input_tasks
    msvc_libgen_task.outputs = [bib_lib_output_file_node]

    Logs.info(f"{msvc_libgen_task.outputs=}")

    return msvc_libgen_task


def run_mvsc_lib_gen(self, task_obj: LibTask):
    Logs.info("Running MSVC lib generation")
    clib_task = task_obj.asterbibc.libtask
    cxxlib_task = task_obj.asterbibcxx.libtask
    fclib_task = task_obj.asterbibfor.libtask
    aster_task = task_obj.asterlib.libtask

    Logs.info(f"Before removal: {clib_task.outputs=}")

    # Lib files are created by MSVC lib generation, so will remove these from the shlib outputs
    clib_task.outputs = [o for o in clib_task.outputs if o.suffix() != ".lib"]
    cxxlib_task.outputs = [o for o in cxxlib_task.outputs if o.suffix() != ".lib"]
    fclib_task.outputs = [o for o in fclib_task.outputs if o.suffix() != ".lib"]
    aster_task.outputs = [o for o in aster_task.outputs if o.suffix() != ".lib"]

    Logs.info(f"After removal: {clib_task.outputs=}")
    Logs.info(f"After removal: {fclib_task.outputs=}")
    Logs.info(f"After removal: {cxxlib_task.outputs=}")
    Logs.info(f"After removal: {aster_task.outputs=}")

    c_input_tasks = [ctask.outputs[0] for ctask in task_obj.asterbibc.tasks]
    cxx_input_tasks = [cxxtask.outputs[0] for cxxtask in task_obj.asterbibcxx.tasks]
    fc_input_tasks = [fctask.outputs[0] for fctask in task_obj.asterbibfor.tasks]
    aster_input_tasks = [ctask.outputs[0] for ctask in task_obj.asterlib.tasks if ctask.outputs[0].suffix() == ".o"]
    Logs.info(f"{aster_input_tasks=}")

    if len(aster_input_tasks) == 0:
        raise Errors.WafError("Failed MSVC lib generation: No aster input tasks found")

    clib_lib_task = create_msvclibgen_task(self, "bibc", c_input_tasks)
    bibcxx_lib_task = create_msvclibgen_task(self, "bibcxx", cxx_input_tasks)
    bibfor_lib_task = create_msvclibgen_task(self, "bibfor", fc_input_tasks)
    bibaster_lib_task = create_msvclibgen_task(self, "aster", aster_input_tasks)

    Logs.info(f"{clib_lib_task.outputs=}")
    Logs.info(f"{fclib_task.outputs=}")
    Logs.info(f"{bibaster_lib_task.outputs=}")

    # filter out all non-lib files
    clib_task_outputs = [x for x in clib_lib_task.outputs if x.suffix() == ".lib"]
    fclib_task_outputs = [x for x in bibfor_lib_task.outputs if x.suffix() == ".lib"]
    bibaster_task_outputs = [x for x in bibaster_lib_task.outputs if x.suffix() == ".lib"]
    bibcxx_task_outputs = [x for x in bibcxx_lib_task.outputs if x.suffix() == ".lib"]

    Logs.info(f"{clib_task_outputs=}")
    Logs.info(f"{fclib_task_outputs=}")
    Logs.info(f"{bibaster_task_outputs=}")
    Logs.info(f"{bibcxx_task_outputs=}")

    clib_task.inputs += bibcxx_task_outputs + fclib_task_outputs
    fclib_task.inputs += bibcxx_task_outputs + clib_task_outputs
    cxxlib_task.inputs += bibaster_task_outputs + clib_task_outputs + fclib_task_outputs
    aster_task.inputs += fclib_task_outputs + clib_task_outputs + bibcxx_task_outputs

    bibc_dll = [o for o in clib_task.outputs if o.suffix() == ".dll"][0]
    bibcxx_dll = [o for o in cxxlib_task.outputs if o.suffix() == ".dll"][0]
    Logs.info(f"{type(bibc_dll)=}{bibc_dll=}")
    Logs.info(f"{type(bibcxx_dll)=}{bibcxx_dll=}")

    mods = ["aster", "aster_core", "aster_fonctions", "med_aster", "libaster"]
    symlink_map = {
        "libaster": bibcxx_dll,
        "aster": bibc_dll,
        "aster_core": bibc_dll,
        "aster_fonctions": bibc_dll,
        "med_aster": bibc_dll,
    }

    aster_lib_dir = pathlib.Path(self.env.ASTERLIBDIR).resolve().absolute()
    input_nodes = []
    output_nodes = []
    for submodule in mods:
        dll_src_node = symlink_map.get(submodule)
        dest = submodule + ".pyd"
        pyd_node = self.bld.bldnode.make_node(dest)
        input_nodes.append(dll_src_node)
        output_nodes.append(pyd_node)
        Logs.info(f"Created symlink: {dll_src_node} -> {pyd_node}")

    msvc_sym_task = self.create_task("msvc_symlink_installer")
    msvc_sym_task.inputs = input_nodes
    msvc_sym_task.outputs = output_nodes
    msvc_sym_task.env = self.env

    Logs.info("Successfully ran MSVC lib generation")


_lib_task_obj: LibTask | None = None
_task_done = False
_compiler_map = {"asterlib": "cxx", "asterbibc": "c", "asterbibfor": "fc", "asterbibcxx": "cxx"}


@TaskGen.feature("cxxshlib", "fcshlib", "cshlib")
@TaskGen.after_method("apply_link", "propagate_uselib_vars")
def set_flags(self) -> None:
    if platform.system() != "Windows":
        return
    name = self.get_name()
    if name == "asterbibc":
        archive_name = "bibc"
    elif name == "asterbibcxx":
        archive_name = "bibcxx"
    elif name == "asterbibfor":
        archive_name = "bibfor"
    elif name == "asterlib":
        archive_name = "aster"
    else:
        Logs.info(f"Skipping {name=}")
        return None
    bld_path = pathlib.Path(self.bld.bldnode.abspath()).resolve().absolute()
    archive_dir = bld_path / archive_name
    args = [f"/LIBPATH:{archive_dir.as_posix()}", f"/WHOLEARCHIVE:{archive_name}.lib", f"{archive_name}.exp"]
    Logs.info(f"Setting flags {args} for {name=}")
    self.link_task.env.append_unique("LINKFLAGS", args)


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
        lib_task_obj.asterbibcxx.tasks + lib_task_obj.asterlib.tasks,
        lib_task_obj.asterbibfor.tasks,
    )
    _task_done = True

