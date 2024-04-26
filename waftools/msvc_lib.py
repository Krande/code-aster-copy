import pathlib
import shutil
from dataclasses import dataclass

from waflib import Logs, TaskGen, Task, Errors, ConfigSet


class msvclibgen(Task.Task):
    run_str = "LIB.exe /OUT:${TGT} ${SRC}"
    color = "BLUE"
    before = ["cshlib", "cxxshlib", "fcshlib"]

    def exec_command(self, cmd, **kw):
        """Execute the command"""
        output_fp = pathlib.Path(self.outputs[0].abspath())
        output_fp.parent.mkdir(parents=True, exist_ok=True)

        clean_name = output_fp.stem.replace("_gen", "")
        root_dir = pathlib.Path(self.generator.bld.root.abspath()).resolve().absolute()
        conda_dir = root_dir / "conda"

        task_gen: TaskGen.task_gen = self.generator
        env: ConfigSet.ConfigSet = task_gen.env
        lib_dir = pathlib.Path(env.LIBDIR)
        root_dir = lib_dir.parent.parent
        libs_dir = root_dir / "libs"

        opts = ["/NOLOGO", "/MACHINE:X64", f"/LIBPATH:{libs_dir}"]

        if clean_name == "bibc":
            bibc_def = conda_dir / "bibc.def"
            opts += [f"/DEF:{bibc_def}"]
        elif clean_name == "bibfor":
            bibfor_def = conda_dir / "bibfor.def"
            opts += [f"/DEF:{bibfor_def}"]

        cmd = cmd[:2] + opts + cmd[2:]

        ret = super().exec_command(cmd, **kw)
        shutil.copy(output_fp, (output_fp.parent / clean_name).with_suffix(".lib"))
        return ret


def build_clang_compilation_db(task_gen, c_tasks, cxx_tasks, fc_tasks=None):
    bld = task_gen.bld
    clang_compilation_database_tasks = c_tasks + cxx_tasks
    database_file = bld.bldnode.make_node("compile_commands.json")
    Logs.info("Build commands will be stored in %s", database_file.path_from(bld.path))
    try:
        root = database_file.read_json()
    except IOError:
        root = []

    bibc_incl = bld.env.INCLUDES + bld.env.INCLUDES_BIBC
    bibcxx_incl = bld.env.INCLUDES + bld.env.INCLUDES_BIBCXX

    Logs.info(f"{bibc_incl=}")
    Logs.info(f"{bibcxx_incl=}")

    args1_cxx = [
        shutil.which("clang-cl.exe"),
        "/nologo",
        "/std:c++17",
        "/FS",
        "/EHsc",
        "/EHc",
        "/permissive-",
        "/Zi",
        "/Od",
        *[f"/I{x}" for x in bibcxx_incl],
        "/DH5_BUILT_AS_DYNAMIC_LIB",
        "/DASTER_HAVE_INTEL_IFORT",
    ]
    args1_c = [
        shutil.which("clang-cl.exe"),
        "/FS",
        "/MD",
        "/nologo",
        "/Zi",
        "/Od",
        *[f"/I{x}" for x in bibc_incl],
        "/DH5_BUILT_AS_DYNAMIC_LIB",
        "/DASTER_HAVE_INTEL_IFORT",
    ]
    top_dir = pathlib.Path(bld.top_dir)
    Logs.info(f"{top_dir=}")
    build_dir = top_dir / f"build/std/{bld.variant}"

    args2 = [
        "/FC",
        "/c",
    ]
    clang_db = dict((x["file"], x) for x in root)
    for task in clang_compilation_database_tasks:
        f_node = task.inputs[0]
        filename = f_node.path_from(task.get_cwd())
        filename_rel = str(filename).replace("..\\..\\", "")
        # Logs.info(f"{filename_rel=}")
        fp = (build_dir / filename).resolve().absolute()
        # Logs.info(f"{fp=}")
        exist_suff = fp.suffix
        dest_o = (build_dir / filename_rel).with_suffix(f"{exist_suff}.1.o")
        # Logs.info(f"{dest_o=}")
        if fp.suffix == ".c":
            args1 = args1_c
        else:
            args1 = args1_cxx
        cmd = args1 + [fp.as_posix()] + args2 + [f"/FoD:{dest_o.as_posix()}"]
        # cmd = task.cmd
        entry = {
            "directory": task.get_cwd().abspath(),
            "arguments": cmd,
            "file": filename,
        }
        clang_db[filename] = entry
    root = list(clang_db.values())
    database_file.write_json(root)
    Logs.info("Compilation database written to %s", database_file.path_from(bld.path))


@dataclass
class TaskObject:
    libtask: Task
    program_task: Task
    tasks: list[Task]
    task_gen: TaskGen


@dataclass
class LibTask:
    asterbibc: TaskObject
    asterpre: TaskObject
    asterbibcxx: TaskObject
    asterbibfor: TaskObject

    def all_tasks_ready(self) -> bool:
        for task in self.__dict__.values():
            if len(task.tasks) == 0:
                Logs.error(f"No tasks found for {task.task_gen.get_name()=}")
                return False
        return True

    def get_missing_tasks(self) -> list[str]:
        missing_tasks = []
        for task in self.__dict__.values():
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
    cxx_aster_object = get_task_object(bld, "asterpre", "cxx")
    cxx_task_object = get_task_object(bld, "asterbibcxx", "cxx")
    fc_task_object = get_task_object(bld, "asterbibfor", "fc")

    return LibTask(c_task_object, cxx_aster_object, cxx_task_object, fc_task_object)


def create_msvclibgen_task(self, lib_name: str, input_tasks) -> Task:
    # Create a task for MSVC lib generation for C
    # the task takes in the outputs of all C tasks
    # it's outputs are bibc.lib and bibc.exp located in the build directory
    bld_path = pathlib.Path(self.bld.bldnode.abspath()).resolve().absolute()
    lib_output_file_path = bld_path / lib_name / f"{lib_name}_gen.lib"

    Logs.info(f"{lib_output_file_path=}")

    # create nodes for the output files
    bib_lib_output_file_node = self.bld.bldnode.make_node(lib_output_file_path.relative_to(bld_path).as_posix())

    msvc_libgen_task = self.create_task("msvclibgen")
    msvc_libgen_task.inputs = input_tasks
    msvc_libgen_task.env = self.env
    msvc_libgen_task.MSVC_LIBGEN_LIB_PATH = lib_output_file_path.as_posix()
    msvc_libgen_task.dep_nodes = input_tasks
    msvc_libgen_task.outputs = [bib_lib_output_file_node]

    return msvc_libgen_task


def run_mvsc_lib_gen(self, task_obj: LibTask):
    Logs.info("Running MSVC lib generation")
    clib_task = task_obj.asterbibc.libtask
    cxxlib_task = task_obj.asterbibcxx.libtask
    fclib_task = task_obj.asterbibfor.libtask

    c_input_tasks = [ctask.outputs[0] for ctask in task_obj.asterbibc.tasks]
    cxx_input_tasks = [cxxtask.outputs[0] for cxxtask in task_obj.asterbibcxx.tasks]
    fc_input_tasks = [fctask.outputs[0] for fctask in task_obj.asterbibfor.tasks]
    aster_input_tasks = [ctask.outputs[0] for ctask in task_obj.asterpre.tasks if ctask.outputs[0].suffix() == ".o"]
    Logs.info(f"{aster_input_tasks=}")

    if len(aster_input_tasks) == 0:
        raise Errors.WafError("Failed MSVC lib generation: No aster input tasks found")

    create_msvclibgen_task(self, "bibc", c_input_tasks)
    bibcxx_lib_task = create_msvclibgen_task(self, "bibcxx", cxx_input_tasks)
    bibfor_lib_task = create_msvclibgen_task(self, "bibfor", fc_input_tasks)
    bibaster_lib_task = create_msvclibgen_task(self, "astertmp", aster_input_tasks)

    clib_task.inputs += bibfor_lib_task.outputs + bibcxx_lib_task.outputs

    Logs.info(f"{clib_task.outputs=}")
    Logs.info(f"{fclib_task.outputs=}")
    Logs.info(f"{bibaster_lib_task.outputs=}")

    # filter out all non-lib files
    clib_task_outputs = [x for x in clib_task.outputs if x.suffix() == ".lib"]
    fclib_task_outputs = [x for x in fclib_task.outputs if x.suffix() == ".lib"]
    bibaster_task_outputs = [x for x in bibaster_lib_task.outputs if x.suffix() == ".lib"]

    Logs.info(f"{clib_task_outputs=}")
    Logs.info(f"{fclib_task_outputs=}")
    Logs.info(f"{bibaster_task_outputs=}")

    fclib_task.inputs += bibcxx_lib_task.outputs + clib_task_outputs
    cxxlib_task.inputs += clib_task_outputs + fclib_task_outputs + bibaster_task_outputs

    Logs.info("Successfully ran MSVC lib generation")


_task_obj = None
_task_done = False
_compiler_map = {"asterpre": "cxx", "asterbibc": "c", "asterbibfor": "fc", "asterbibcxx": "cxx"}


@TaskGen.feature("cxxshlib", "fcshlib", "cshlib", "cxx")
@TaskGen.after_method("apply_link")
def make_msvc_modifications(self: TaskGen.task_gen):
    name = self.get_name()
    Logs.info(f"task: {name=}")

    if name not in _compiler_map.keys():
        return

    global _task_obj
    global _task_done
    if _task_done:
        return

    if _task_obj is None:
        task_obj = extract_main_tasks(self)
        _task_obj = task_obj
    else:
        task_obj = _task_obj
        if not task_obj.all_tasks_ready():
            missing_task_names = task_obj.get_missing_tasks()
            Logs.error(f"Missing tasks before: {missing_task_names}")
            aster_object = get_task_object(self.bld, name, _compiler_map[name])
            if name == "asterpre":
                aster_object.tasks = self.tasks
            setattr(task_obj, name, aster_object)

    if not task_obj.all_tasks_ready():
        missing_task_names = task_obj.get_missing_tasks()
        Logs.error(f"Missing tasks @{name} after: {missing_task_names}")
    else:
        Logs.info("All tasks are ready")
        run_mvsc_lib_gen(self, task_obj)

        build_clang_compilation_db(
            self,
            task_obj.asterbibc.tasks,
            task_obj.asterbibcxx.tasks,
            task_obj.asterbibfor.tasks,
        )
        _task_done = True
