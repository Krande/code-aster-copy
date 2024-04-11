import os
import pathlib
import shutil
from dataclasses import dataclass

from waflib import Logs, TaskGen, Task, Node, Utils, Errors


class msvclibgen(Task.Task):
    run_str = "lib /OUT:{MSVC_LIBGEN_LIB_PATH} @{MSVC_LIBGEN_INPUT_FILE_TXT}"
    color = "BLUE"

    def run(self):
        """Execute the command"""
        kwargs = self.env.get_merged_dict()
        kwargs["cwd"] = self.get_cwd().abspath()
        Logs.info()

        output_fp = pathlib.Path(self.outputs[0].abspath())
        output_fp.parent.mkdir()
        ofiles = output_fp.parent / f"{output_fp.stem}_input_files.txt"
        with open(ofiles, 'w') as f:
            for i in self.inputs:
                f.write(f"{i.abspath()}\n")
        cmds = self.run_str.format(MSVC_LIBGEN_LIB_PATH=output_fp.as_posix(),
                                   MSVC_LIBGEN_INPUT_FILE_TXT=ofiles.as_posix())

        ret, stdout, stderr = Utils.run_process(cmds, kwargs)
        if ret:
            error = Errors.WafError("Command %r returned %r" % (cmds, ret))
            error.returncode = ret
            error.stderr = stderr
            error.stdout = stdout
            raise error

    def exec_command(self, cmd, **kw):
        return super().exec_command(cmd, **kw)


def scrape_all_tasks(self):
    """Scrape all outputs from all tasks."""
    all_outputs = []
    lib_deps = {}
    for mod in ["asterbibfor", "asterbibc", "asterbibcxx"]:
        for task in self.bld.get_tgen_by_name(mod).tasks:
            if task.__class__.__name__ not in ("c", "cxx", "fc"):
                continue
            fp_in = pathlib.Path(task.inputs[0].abspath())
            lib_deps[fp_in] = task
            all_outputs.extend(task.outputs)
    return all_outputs, lib_deps


def fix_specific_c_bibfor_includes(self):
    config_dir = pathlib.Path(__file__).parent.parent / "config"
    cshlib_txt = config_dir / "cshlib.txt"
    cxxshlib_txt = config_dir / "cxxshlib.txt"

    added_links = {
        "cshlib": cshlib_txt.read_text().splitlines(),
        "cxxshlib": cxxshlib_txt.read_text().splitlines(),
    }
    top_dir = pathlib.Path(self.bld.top_dir)
    all_outputs, lib_deps = scrape_all_tasks(self)
    # Logs.info(f"{all_outputs=}, {lib_deps=}")
    # Logs.info(f"{self.bld.env.INCLUDES_BIBFOR=}")

    # fscan = fc_scan.fortran_parser(self.bld.env.INCLUDES_BIBFOR)
    # fp1 = top_dir / added_links["cshlib"][1]
    # dep_task = lib_deps[fp1]
    # n = self.path.find_node(dep_task.inputs[0])
    # Logs.info(f"{n}")
    # fscan.start(n)
    # Logs.info(f"{fscan.nodes=}")

    for task in self.bld.get_tgen_by_name("asterbibc").tasks:
        tname = task.__class__.__name__
        if tname != "cshlib":
            continue
        deps = added_links[tname]
        for dep in deps:
            dep_fp = top_dir / dep
            if dep_fp not in lib_deps:
                Logs.error(f"Dependency {dep_fp} not found for bibc")
                continue
            dep_task = lib_deps[dep_fp]
            task.inputs.extend(dep_task.outputs)
            Logs.info(f"Adding {dep_task} to {task}")
            # Logs.info(f"{task.inputs=}")

    for task in self.bld.get_tgen_by_name("asterbibcxx").tasks:
        tname = task.__class__.__name__
        if tname != "cxxshlib":
            continue
        # task.inputs.extend(all_outputs)
        deps = added_links[tname]
        for dep in deps:
            dep_fp = top_dir / dep
            if dep_fp not in lib_deps:
                Logs.error(f"Dependency {dep_fp} not found for bibcxx")
                continue
            dep_task = lib_deps[dep_fp]
            task.inputs.extend(dep_task.outputs)
            Logs.info(f"Adding {dep_task} to {task}")
            # Logs.info(f"{task.inputs=}")


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
class LibTask:
    clib_task: Task
    cxxlib_task: Task
    fclib_task: Task

    c_program_task: Task
    cxx_program_task: Task
    fc_program_task: Task

    c_tasks: list[Task]
    cxx_tasks: list[Task]
    fc_tasks: list[Task]


def extract_main_tasks(self) -> LibTask:
    fc_task = None
    fc_program_task = None
    fc_tasks = []
    for task in self.bld.get_tgen_by_name("asterbibfor").tasks:
        if task.__class__.__name__ == "fcshlib":
            fc_task = task
        if task.__class__.__name__ == "fcprogram":
            fc_program_task = task
        if task.__class__.__name__ == "fc":
            fc_tasks.append(task)
        # Logs.info(f"{task.__class__.__name__=}")
        # Logs.info(f"{task.outputs=}")
        # Logs.info(f"{task.inputs=}")
        # Logs.info(f"{task=}")

    c_task = None
    c_program_task = None
    c_tasks = []
    for task in self.bld.get_tgen_by_name("asterbibc").tasks:
        if task.__class__.__name__ == "cshlib":
            c_task = task
        if task.__class__.__name__ == "cprogram":
            c_program_task = task
        if task.__class__.__name__ == "c":
            c_tasks.append(task)
        # Logs.info(f"{task.__class__.__name__=}")
        # Logs.info(f"{task.outputs=}")
        # Logs.info(f"{task.inputs=}")
        # Logs.info(f"{task=}")

    cxx_task = None
    cxx_program_task = None
    cxx_tasks = []
    for task in self.bld.get_tgen_by_name("asterbibcxx").tasks:
        if task.__class__.__name__ == "cxxshlib":
            cxx_task = task
        if task.__class__.__name__ == "cxxprogram":
            cxx_program_task = task
        if task.__class__.__name__ == "cxx":
            cxx_tasks.append(task)
        # Logs.info(f"{task.__class__.__name__=}")
        # Logs.info(f"{task.outputs=}")
        # Logs.info(f"{task.inputs=}")
        # Logs.info(f"{task=}")

    return LibTask(c_task, cxx_task, fc_task, c_program_task, cxx_program_task, fc_program_task, c_tasks, cxx_tasks,
                   fc_tasks)


def force_bibfor_seq(task_obj: LibTask):
    Logs.info("Forcing linking order: bibfor -> bibc -> bibcxx")

    c_task = task_obj.clib_task
    cxx_task = task_obj.cxxlib_task
    fc_task = task_obj.fclib_task

    c_task.dep_nodes.extend(fc_task.outputs)
    cxx_task.dep_nodes.extend(fc_task.outputs)
    cxx_task.dep_nodes.extend(c_task.outputs)

    extra_c_deps = {
        "envima.c",
        "aster_module.c",
        "aster_core_module.c",
        "aster_utils.c",
        "shared_vars.c",
        "aster_exceptions.c",
        "utflsh.c",
        "uttrst.c",
        "iodr.c",
        "rmfile.c",
        "hpalloc.c",
        "inisig.c",
        "aster_mpi.c",
        "indik8.c",
        "matfpe.c",
        "erfcam.c",
        "hpdeallc.c",
        "dll_umat.c",
        "hdftsd.c",
        "dll_interface.c",
        "dll_register.c",
        "uttcsm.c",
        "fetsco.c",
        "gpmetis_aster.c",
        "libinfos.c",
        "kloklo.c",
        "cpfile.c",
        "hdfrsv.c",
        "datetoi.c",
        "strmov.c",
        "hanfpe.c",
        "mempid.c"
    }
    found = set()
    # need to add output of som "*.c" files to bibfor to resolve shared c symbols
    for task in task_obj.c_tasks:
        abs_path = pathlib.Path(task.inputs[0].abspath()).resolve()
        for extra_dep in extra_c_deps:
            if extra_dep == abs_path.name:
                fc_task.inputs.extend(task.outputs)
                found.add(extra_dep)
                break

    # check if all extra deps are found
    if found != extra_c_deps:
        Logs.error(f"Extra deps not found: {extra_c_deps - found}")

    # Not quite sure about these lines
    # c_task.inputs.extend(fc_task.outputs)
    # cxx_task.inputs.extend(fc_task.outputs)
    # cxx_task.inputs.extend(c_task.outputs)

    # Logs.info(f"{c_task.priority()=}")
    # Logs.info(f"{cxx_task.priority()=}")
    # Logs.info(f"{fc_task.priority()=}")


def create_msvclibgen_task(self, lib_name: str, input_tasks) -> Task:
    # Create a task for MSVC lib generation for C
    # the task takes in the outputs of all C tasks
    # it's outputs are bibc.lib and bibc.exp located in the build directory
    bld_path = pathlib.Path(self.bld.bldnode.abspath()).resolve().absolute()
    bibc_lib_output_file_path = bld_path / lib_name / f"{lib_name}.lib"
    Logs.info(f"{bibc_lib_output_file_path=}")

    # create nodes for the output files
    bib_lib_output_file_node = self.bld.bldnode.make_node(bibc_lib_output_file_path.relative_to(bld_path).as_posix())

    msvc_libgen_task = self.create_task("msvclibgen")
    msvc_libgen_task.inputs = input_tasks
    msvc_libgen_task.env = self.env
    msvc_libgen_task.MSVC_LIBGEN_LIB_PATH = bibc_lib_output_file_path.as_posix()
    msvc_libgen_task.dep_nodes = input_tasks
    msvc_libgen_task.outputs = [bib_lib_output_file_node]

    return msvc_libgen_task


def run_mvsc_lib_gen(self, task_obj: LibTask):
    Logs.info("Running MSVC lib generation")
    clib_task = task_obj.clib_task
    cxxlib_task = task_obj.cxxlib_task
    fclib_task = task_obj.fclib_task

    # c_input_tasks = [ctask.outputs[0] for ctask in task_obj.c_tasks]
    cxx_input_tasks = [cxxtask.outputs[0] for cxxtask in task_obj.cxx_tasks]
    fc_input_tasks = [fctask.outputs[0] for fctask in task_obj.fc_tasks]

    # bibc_output = create_msvclibgen_task(self, "bibc", c_input_tasks)
    bibcxx_lib_task = create_msvclibgen_task(self, "bibcxx", cxx_input_tasks)
    bibfor_lib_task = create_msvclibgen_task(self, "bibfor", fc_input_tasks)

    clib_task.inputs += bibfor_lib_task.outputs + bibcxx_lib_task.outputs
    # clib_task.run_str = ""

    # cxxlib_task.inputs = bibc_output + bibfor_output
    # cxxlib_task.run_str = ""

    # fclib_task.inputs = bibc_output + bibcxx_output
    # fclib_task.run_str = "${LINK_FC} ${FCLINKFLAGS} ${LINKFLAGS} ${FCLNK_SRC_F}${SRC} ${FCLNK_TGT_F}${TGT[0].abspath()} ${RPATH_ST:RPATH} ${FCSTLIB_MARKER} ${FCSTLIBPATH_ST:STLIBPATH} ${FCSTLIB_ST:STLIB} ${FCSHLIB_MARKER} ${FCLIBPATH_ST:LIBPATH} ${FCLIB_ST:LIB} ${LDFLAGS}"
    # cxxlib_task.inputs += clib_task.outputs
    # cxxlib_task.outputs = [cxxlib_task.outputs[0].change_ext(".lib")]
    # fclib_task.outputs = [fclib_task.outputs[0].change_ext(".lib")]




@TaskGen.feature("cshlib")
@TaskGen.after_method("apply_link")
def make_msvc_modifications(self):
    task_obj = extract_main_tasks(self)

    # Print depends order for each task
    # Logs.info(f"{task_obj.fc_task.dep_nodes=}")
    # Logs.info(f"{task_obj.c_task.dep_nodes=}")
    # Logs.info(f"{task_obj.cxx_task.dep_nodes=}")

    # These will Force linking order: bibfor -> bibc -> bibcxx
    if os.getenv("FORCE_BIBFOR_SEQUENCE"):
        force_bibfor_seq(task_obj)

    if os.getenv("MANUALLY_ADD_BIBFOR_DEPS"):
        fix_specific_c_bibfor_includes(self)

    run_mvsc_lib_gen(self, task_obj)

    Logs.info(f"{task_obj.clib_task.priority()=}")
    Logs.info(f"{task_obj.cxxlib_task.priority()=}")
    Logs.info(f"{task_obj.fclib_task.priority()=}")

    build_clang_compilation_db(self, task_obj.c_tasks, task_obj.cxx_tasks)
