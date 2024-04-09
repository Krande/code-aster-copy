import os
import pathlib
import shutil

from waflib import Logs, TaskGen


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
    build_dir = top_dir / "build/std/debug"

    args2 = [
        "/FC",
        "/c",
    ]
    clang_db = dict((x["file"], x) for x in root)
    for task in clang_compilation_database_tasks:
        f_node = task.inputs[0]
        filename = f_node.path_from(task.get_cwd())
        filename_rel = str(filename).replace("..\\..\\", '')
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


@TaskGen.feature("cshlib")
@TaskGen.after_method("apply_link")
def def_prep_fc_c_linking(self):
    fc_task = None
    for task in self.bld.get_tgen_by_name("asterbibfor").tasks:
        if task.__class__.__name__ != "fcshlib":
            continue
        Logs.info(f"{task.__class__.__name__=}")
        Logs.info(f"{task.outputs=}")
        # Logs.info(f"{task.inputs=}")
        # Logs.info(f"{task=}")
        fc_task = task

    c_task = None
    c_tasks = []
    for task in self.bld.get_tgen_by_name("asterbibc").tasks:
        if task.__class__.__name__ != "cshlib":
            if task.__class__.__name__ == "c":
                c_tasks.append(task)
            continue
        Logs.info(f"{task.__class__.__name__=}")
        Logs.info(f"{task.outputs=}")
        # Logs.info(f"{task.inputs=}")
        # Logs.info(f"{task=}")
        c_task = task

    cxx_task = None
    cxx_tasks = []
    for task in self.bld.get_tgen_by_name("asterbibcxx").tasks:
        if task.__class__.__name__ != "cxxshlib":
            if task.__class__.__name__ == "cxx":
                cxx_tasks.append(task)
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
    if os.getenv("FORCE_BIBFOR_SEQUENCE"):
        Logs.info("Forcing linking order: bibfor -> bibc -> bibcxx")
        c_task.dep_nodes.extend(fc_task.outputs)
        cxx_task.dep_nodes.extend(fc_task.outputs)
        cxx_task.dep_nodes.extend(c_task.outputs)
        envima_c_task = None
        # need to add output of "envima.c" to bibfor to resolve shared c symbols
        for task in c_tasks:
            if "envima.c" in task.inputs[0].abspath():
                envima_c_task = task
                break
        if envima_c_task:
            fc_task.inputs.extend(envima_c_task.outputs)
        else:
            Logs.error("envima.c output not found in c_tasks")

        # Not quite sure about these lines
        # c_task.inputs.extend(fc_task.outputs)
        # cxx_task.inputs.extend(fc_task.outputs)
        # cxx_task.inputs.extend(c_task.outputs)

    if os.getenv("MANUALLY_ADD_BIBFOR_DEPS"):
        fix_specific_c_bibfor_includes(self)

    Logs.info(f"{c_task.priority()=}")
    Logs.info(f"{cxx_task.priority()=}")
    Logs.info(f"{fc_task.priority()=}")

    build_clang_compilation_db(self, c_tasks, cxx_tasks)
