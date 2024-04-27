import pathlib
import shutil

from waflib import Logs


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
