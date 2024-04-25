import os
import pathlib
import shutil
from dataclasses import dataclass

from waflib import Logs, TaskGen, Task, Node, Utils, Errors


class msvclibgen(Task.Task):
    run_str = "${LIB_EXE} /OUT:${MSVC_LIBGEN_LIB_PATH} @${MSVC_LIBGEN_INPUT_FILE_TXT}"
    color = "BLUE"

    def exec_command(self, cmd, **kw):
        """Execute the command"""
        environ = None
        if "env" in kw:
            environ = kw["env"]
            del kw["env"]

        if environ is None:
            environ = os.environ.copy()

        output_fp = pathlib.Path(self.outputs[0].abspath())

        ofiles = output_fp.parent / f"{output_fp.stem}_input_files.txt"
        ofiles.parent.mkdir(parents=True, exist_ok=True)
        with open(ofiles, "w") as f:
            for i in self.inputs:
                f.write(f"{i.abspath()}\n")

        environ["MSVC_LIBGEN_LIB_PATH"] = output_fp.as_posix()
        environ["MSVC_LIBGEN_INPUT_FILE_TXT"] = ofiles.as_posix()

        kw["env"] = environ
        kw["shell"] = True
        # Do this for now.
        cmd = f"LIB.exe /OUT:{output_fp.as_posix()} @{ofiles.as_posix()}".split()
        Logs.info("Executing %r", cmd)
        clean_name = output_fp.stem.replace("_gen", "")
        ret = super().exec_command(cmd, **kw)
        shutil.copy(output_fp, (output_fp.parent / clean_name).with_suffix(".lib"))
        return ret


def build_clang_compilation_db(task_gen, c_tasks, cxx_tasks):
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


@dataclass
class LibTask:
    c_task: TaskObject
    c_aster_task: TaskObject
    cxx_task: TaskObject
    fc_task: TaskObject


def get_task_object(bld: TaskGen, taskgen_name: str, compiler_prefix: str) -> TaskObject:
    lib_task = None
    program_task = None
    tasks = []
    for task in bld.get_tgen_by_name(taskgen_name).tasks:
        if task.__class__.__name__ == f"{compiler_prefix}shlib":
            lib_task = task
        if task.__class__.__name__ == f"{compiler_prefix}program":
            program_task = task
        if task.__class__.__name__ == compiler_prefix:
            tasks.append(task)

    if len(tasks) == 0:
        Logs.info(f"No tasks found for {taskgen_name}")

    return TaskObject(lib_task, program_task, tasks)


def extract_main_tasks(bld) -> LibTask:
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
    clib_task = task_obj.c_task.libtask
    cxxlib_task = task_obj.cxx_task.libtask
    fclib_task = task_obj.fc_task.libtask

    cxx_input_tasks = [cxxtask.outputs[0] for cxxtask in task_obj.cxx_task.tasks]
    fc_input_tasks = [fctask.outputs[0] for fctask in task_obj.fc_task.tasks]

    bibcxx_lib_task = create_msvclibgen_task(self, "bibcxx", cxx_input_tasks)
    bibfor_lib_task = create_msvclibgen_task(self, "bibfor", fc_input_tasks)

    clib_task.inputs += bibfor_lib_task.outputs + bibcxx_lib_task.outputs

    Logs.info(f"{clib_task.outputs=}")
    Logs.info(f"{fclib_task.outputs=}")

    # filter out all non-lib files
    clib_task_outputs = [x for x in clib_task.outputs if x.suffix() == ".lib"]
    fclib_task_outputs = [x for x in fclib_task.outputs if x.suffix() == ".lib"]

    Logs.info(f"{clib_task_outputs=}")
    Logs.info(f"{fclib_task_outputs=}")

    fclib_task.inputs += bibcxx_lib_task.outputs + clib_task_outputs
    cxxlib_task.inputs += clib_task_outputs + fclib_task_outputs

    bld_path = pathlib.Path(self.bld.bldnode.abspath()).resolve().absolute()
    lib_output_file_path = bld_path / "aster" / f"py_aster_gen.lib"
    bib_lib_output_file_node = self.bld.bldnode.make_node(lib_output_file_path.relative_to(bld_path).as_posix())

    py_c_obj = self.path.find_resource("supervis/python.c.2.o")
    input_tasks = [py_c_obj]

    msvc_libgen_task = self.create_task("msvclibgen")
    msvc_libgen_task.inputs = input_tasks
    msvc_libgen_task.env = self.env
    msvc_libgen_task.MSVC_LIBGEN_LIB_PATH = lib_output_file_path.as_posix()
    msvc_libgen_task.dep_nodes = input_tasks
    msvc_libgen_task.outputs = [bib_lib_output_file_node]


@TaskGen.feature("cshlib")
@TaskGen.after_method("apply_link")
def make_msvc_modifications_pre(self):
    task_obj = extract_main_tasks(self.bld)
    py_c = self.path.find_resource("supervis/python.c")
    self.bld.add_manual_dependency(
        py_c,
        self.path.find_resource("supervis/aster_core_module.c")
    )

    c_input_tasks = [ctask.outputs[0] for ctask in task_obj.c_task.tasks]
    create_msvclibgen_task(self, "bibc", c_input_tasks)


@TaskGen.feature("cshlib")
@TaskGen.after_method("apply_link")
def make_msvc_modifications(self):
    task_obj = extract_main_tasks(self.bld)

    run_mvsc_lib_gen(self, task_obj)

    Logs.info(f"{task_obj.c_task.libtask.priority()=}")
    Logs.info(f"{task_obj.cxx_task.libtask.priority()=}")
    Logs.info(f"{task_obj.fc_task.libtask.priority()=}")

    build_clang_compilation_db(self, task_obj.c_task.tasks, task_obj.cxx_task.tasks)
