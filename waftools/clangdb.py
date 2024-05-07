import pathlib
import shutil

from waflib import Logs, Task


def build_clang_compilation_db(task_gen, c_tasks, cxx_tasks, fc_tasks=None):
    bld = task_gen.bld
    clang_compilation_database_tasks = c_tasks + cxx_tasks + (fc_tasks or [])
    database_file = bld.bldnode.make_node("compile_commands.json")
    Logs.info("Build commands will be stored in %s", database_file.path_from(bld.path))
    try:
        root = database_file.read_json()
    except IOError:
        root = []

    top_dir = pathlib.Path(bld.top_dir).resolve().absolute()
    build_dir = top_dir / f"build/std/{bld.variant}"

    clang_db = dict((x["file"], x) for x in root)
    for task in clang_compilation_database_tasks:
        task: Task.Task
        f_node = task.inputs[0]
        filename = f_node.path_from(task.get_cwd())
        filename_rel = str(filename).replace("..\\..\\", "")
        fp = (build_dir / filename).resolve().absolute()
        exist_suff = fp.suffix
        dest_o = (build_dir / filename_rel).with_suffix(f"{exist_suff}.1.o")

        args = []
        for define in task.env.DEFINES:
            args.append(f"/D{define}")
        for incl in task.env.INCLUDES:
            args.append(f"/I{incl}")

        if fp.suffix.lower() == ".c":
            compiler = task.env.CC
            args += task.env.CFLAGS
        elif fp.suffix.lower() == ".cxx":
            compiler = task.env.CXX
            args += task.env.CXXFLAGS
        elif fp.suffix.lower() == ".f90":
            compiler = task.env.FC
            args += task.env.FCFLAGS
        else:
            raise ValueError(f"Unknown file extension: {fp.suffix} in {fp}")

        compiler = [pathlib.Path(shutil.which(compiler[0])).resolve().absolute().as_posix()]
        if fp.suffix.lower() == ".f90":
            cmd = compiler + args + ["/o", f"{dest_o.as_posix()}"] + [fp.as_posix()]
        else:
            cmd = compiler + args + [fp.as_posix()] + [f"/Fo{dest_o.as_posix()}"]

        entry = {
            "directory": task.get_cwd().abspath(),
            "arguments": cmd,
            "file": filename,
        }
        clang_db[filename] = entry

    root = list(clang_db.values())
    database_file.write_json(root)
    Logs.info("Compilation database written to %s", database_file.path_from(bld.path))
