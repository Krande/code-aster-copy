import pathlib
import shutil
import time
import os
import sys
import json
import subprocess

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

    total = len(clang_compilation_database_tasks)
    if total:
        Logs.info("Preparing compilation database entries: %d tasks", total)

    # Cache compiler path lookups to avoid repeated shutil.which on large projects
    compiler_cache = {}

    def resolve_compiler(bin_name_list):
        if not bin_name_list:
            return []
        key = bin_name_list[0]
        val = compiler_cache.get(key)
        if val is None:
            resolved = pathlib.Path(shutil.which(key)).resolve().absolute().as_posix()
            val = [resolved]
            compiler_cache[key] = val
        return val

    clang_db = dict((x["file"], x) for x in root)

    start = time.perf_counter()
    next_progress = 0.1  # 10% increments

    for idx, task in enumerate(clang_compilation_database_tasks, 1):
        task: Task.Task
        f_node = task.inputs[0]
        filename = f_node.path_from(task.get_cwd())
        filename_rel = str(filename).replace("..\\..\\", "")
        fp = (build_dir / filename).resolve().absolute()
        exist_suff = fp.suffix
        dest_o = (build_dir / filename_rel).with_suffix(f"{exist_suff}.1.o")

        args = []
        for define in getattr(task.env, "DEFINES", []):
            args.append(f"/D{define}")
        for incl in getattr(task.env, "INCLUDES", []):
            args.append(f"/I{incl}")

        suff = fp.suffix.lower()
        if suff == ".c":
            compiler = resolve_compiler(getattr(task.env, "CC", []))
            args += getattr(task.env, "CFLAGS", [])
        elif suff == ".cxx":
            compiler = resolve_compiler(getattr(task.env, "CXX", []))
            args += getattr(task.env, "CXXFLAGS", [])
        elif suff == ".f90":
            compiler = resolve_compiler(getattr(task.env, "FC", []))
            args += getattr(task.env, "FCFLAGS", [])
        else:
            raise ValueError(f"Unknown file extension: {fp.suffix} in {fp}")

        if suff == ".f90":
            cmd = compiler + args + ["/o", f"{dest_o.as_posix()}"] + [fp.as_posix()]
        else:
            cmd = compiler + args + [fp.as_posix()] + [f"/Fo{dest_o.as_posix()}"]

        entry = {
            "directory": task.get_cwd().abspath(),
            "arguments": cmd,
            "file": filename,
        }
        clang_db[filename] = entry

        # Progress log every 10% (and for very large sets at most every 1000 tasks)
        if total:
            ratio = idx / total
            if ratio >= next_progress or (idx % 1000 == 0):
                Logs.info("Preparing compilation database: %d/%d (%.0f%%)", idx, total, ratio * 100)
                next_progress += 0.1

    elapsed = time.perf_counter() - start
    root = list(clang_db.values())

    # Async option: offload the final JSON write to a background subprocess
    if os.environ.get("ASTER_CLANGDB_MODE", "").lower() == "async":
        snapshot_node = bld.bldnode.make_node("compile_commands.snapshot.json")
        snapshot_path = snapshot_node.abspath()
        out_path = database_file.abspath()
        # Write snapshot atomically
        tmp_snap = snapshot_path + ".tmp"
        with open(tmp_snap, "w", encoding="utf-8") as f:
            json.dump(root, f)
        os.replace(tmp_snap, snapshot_path)

        # Build subprocess command to finalize JSON (atomic write + log)
        cmd = [sys.executable, "-m", "waftools.clangdb", "--from-snapshot", snapshot_path, "--out", out_path]
        creationflags = 0
        popen_kwargs = {}
        if os.name == "nt":
            # Detached background process on Windows
            creationflags = getattr(subprocess, "DETACHED_PROCESS", 0) | getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
            popen_kwargs["creationflags"] = creationflags
        else:
            popen_kwargs["start_new_session"] = True

        try:
            subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, **popen_kwargs)
            Logs.info(
                "Compilation database generation continues in background (async). Snapshot has %d entries; elapsed so far: %.2f s",
                len(root),
                elapsed,
            )
        except Exception as e:
            Logs.warn("Failed to spawn async clangdb writer (%s); falling back to synchronous write", e)
            database_file.write_json(root)
            Logs.info(
                "Compilation database written to %s (%d entries in %.2f s)",
                database_file.path_from(bld.path),
                len(root),
                elapsed,
            )
        return

    # Synchronous path (default)
    database_file.write_json(root)
    Logs.info(
        "Compilation database written to %s (%d entries in %.2f s)",
        database_file.path_from(bld.path),
        len(root),
        elapsed,
    )


if __name__ == "__main__":
    # CLI for async writer: python -m waftools.clangdb --from-snapshot <snapshot> --out <output>
    import argparse
    _ap = argparse.ArgumentParser(description="Finalize compile_commands.json from snapshot.")
    _ap.add_argument("--from-snapshot", dest="snapshot", required=True, help="Path to snapshot JSON")
    _ap.add_argument("--out", dest="out", required=True, help="Path to final compile_commands.json")
    _args = _ap.parse_args()

    t0 = time.perf_counter()
    with open(_args.snapshot, "r", encoding="utf-8") as _f:
        _data = json.load(_f)

    _tmp = _args.out + ".tmp"
    with open(_tmp, "w", encoding="utf-8") as _f:
        json.dump(_data, _f)
    os.replace(_tmp, _args.out)

    _elapsed = time.perf_counter() - t0
    try:
        from waflib import Logs as _Logs  # type: ignore
        _Logs.info(
            "Async compilation database written to %s (%d entries in %.2f s)",
            _args.out,
            len(_data),
            _elapsed,
        )
    except Exception:
        print(
            f"Async compilation database written to {_args.out} ({len(_data)} entries in {_elapsed:.2f} s)",
            file=sys.stdout,
        )
