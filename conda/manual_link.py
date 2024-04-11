import os
import pathlib
import subprocess

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent
THIS_DIR = pathlib.Path(__file__).resolve().parent
BUILD_DIR = ROOT_DIR / "build" / "std" / "debug"

CONDA_PREFIX = os.getenv("CONDA_PREFIX")


def create_args(
        lib_name, dep_files: list[str], extra_deps: list[str], extra_flags: list[str] = None,
        pre_flags: list[str] = None
) -> list[str]:
    extra_flags = extra_flags or []
    pre_flags = pre_flags or []
    return [
        "call_link.bat",
        "/nologo",
        "/MANIFEST",
        "/subsystem:console",
        f"/IMPLIB:{lib_name}\\{lib_name}.lib",
        "/DLL",
        f"/LIBPATH:{CONDA_PREFIX}/libs",
        f"/LIBPATH:{CONDA_PREFIX}/include",
        f"/LIBPATH:{ROOT_DIR}\\{lib_name}\\include",
        f"/LIBPATH:{CONDA_PREFIX}/Library/lib",
        f"/LIBPATH:{CONDA_PREFIX}/Library/bin",
        *pre_flags,
        *dep_files,
        f"/OUT:{BUILD_DIR}\\{lib_name}\\{lib_name}.dll",
        *extra_deps,
        *extra_flags,
    ]


def filter_objects(objs: set[pathlib.Path], filter_func) -> set[pathlib.Path]:
    return set(filter(filter_func, objs))


def bibc(include_all_bibfor=False):
    bibc_objects = set(
        (BUILD_DIR / x).absolute() for x in (THIS_DIR / "bibc_default_order.txt").read_text().splitlines()
    )
    bibfor_objects = set(
        (BUILD_DIR / f"{x}.1.o").absolute() for x in (THIS_DIR / "cshlib.txt").read_text().splitlines()
    )

    if include_all_bibfor:
        all_bibfor_files = (BUILD_DIR / "bibfor").rglob("*.o")
        bibfor_objects.update(set(x.absolute() for x in all_bibfor_files))

    cleaned_bibfor_objects = filter_objects(bibfor_objects, lambda x: not x.name.startswith("ar_d"))
    # print the removed objects
    removed_objects = bibfor_objects - cleaned_bibfor_objects
    print("Removed objects:")
    print("\n".join(map(str, removed_objects)))

    bibc_obj_list_path = ROOT_DIR / "tmp_bibc_objects.txt"

    extra_deps = [
        "esmumps.lib",
        "smumps.lib",
        "scotch.lib",
        "scotcherr.lib",
        "metis.lib",
        "medC.lib",
        "hdf5.lib",
        "pthread.lib",
        "z.lib",
        "bibfor.lib",
        "bibcxx.lib",
    ]
    input_args = create_args("bibc", [f"@{bibc_obj_list_path.name}"], extra_deps, ["/DEBUG", "/INCREMENTAL:NO"])
    input_args += [f"/LIBPATH:{BUILD_DIR}/bibc", f"/LIBPATH:{BUILD_DIR}/bibfor", f"/LIBPATH:{BUILD_DIR}/bibcxx"]

    with open(bibc_obj_list_path, "w") as f:
        f.write("\n".join(map(str, bibc_objects)))
        # f.write("\n")
        # f.write("\n".join(map(str, cleaned_bibfor_objects)))

    # Now call that batch file
    print(" ".join(input_args))
    subprocess.run(input_args, shell=True, cwd=ROOT_DIR)


def bibfor():
    all_bibfor_files = (BUILD_DIR / "bibfor").rglob("*.o")
    bibfor_objects = set(x.absolute() for x in all_bibfor_files)

    cleaned_bibfor_objects = filter_objects(bibfor_objects, lambda x: not x.name.startswith("ar_d"))
    bibfor_obj_list_path = ROOT_DIR / "tmp_bibfor_objects.txt"

    # print the removed objects
    removed_objects = bibfor_objects - cleaned_bibfor_objects
    print("Removed objects:")
    print("\n".join(map(str, removed_objects)))
    cleaned_bibfor_objects_sorted = sorted(cleaned_bibfor_objects, key=lambda x: x.name)

    with open(bibfor_obj_list_path, "w") as f:
        f.write("\n".join(map(str, cleaned_bibfor_objects_sorted)))

    dep_files = [f"@{bibfor_obj_list_path.name}"]
    extra_deps = [
        "dmumps.lib",
        "zmumps.lib",
        "smumps.lib",
        "cmumps.lib",
        "esmumps.lib",
        "mumps_common.lib",
        "pord.lib",
        "pthread.lib",
        "mkl_intel_lp64_dll.lib",
        "mkl_intel_thread_dll.lib",
        "mkl_core_dll.lib"
    ]
    extra_flags = [
        # "/NOENTRY"
        # "/FORCE",
        # "/WHOLEARCHIVE"
    ]

    input_args = create_args("bibfor", dep_files, extra_deps, extra_flags=extra_flags)

    subprocess.run(input_args, shell=True, cwd=ROOT_DIR)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Manually Link Code Aster libraries")
    parser.add_argument("--bibc", action="store_true", help="Link the bibc library")
    parser.add_argument("--bibfor", action="store_true", help="Link the bibfor library")
    parser.add_argument("--dll-to-lib", help="Create a .lib file from a .dll file")

    args = parser.parse_args()

    if args.bibc:
        bibc()
    if args.bibfor:
        bibfor()
    if args.dll_to_lib:
        create_lib_from_dll(args.dll_to_lib)
