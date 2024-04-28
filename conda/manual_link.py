import os
import pathlib

from msvc_utils import call_using_env

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent
THIS_DIR = pathlib.Path(__file__).resolve().parent
BUILD_DIR = ROOT_DIR / "build" / "std" / "debug"
TMP_DIR = THIS_DIR / "temp"
TMP_DIR.mkdir(exist_ok=True)

CONDA_PREFIX = os.getenv("CONDA_PREFIX")


# https://learn.microsoft.com/en-us/cpp/build/reference/linker-options?view=msvc-170
def get_default_flags():
    return [
        "/nologo",
        "/MANIFEST",
        "/subsystem:console",
        "/DLL",
    ]


def get_default_lib_paths(lib_name: str):
    return [
        f"/LIBPATH:{CONDA_PREFIX}/libs",
        f"/LIBPATH:{CONDA_PREFIX}/include",
        f"/LIBPATH:{ROOT_DIR}/{lib_name}/include",
        f"/LIBPATH:{CONDA_PREFIX}/Library/lib",
        f"/LIBPATH:{CONDA_PREFIX}/Library/bin",
        f"/LIBPATH:{TMP_DIR}",
    ]


def filter_objects(objs: set[pathlib.Path], filter_func) -> set[pathlib.Path]:
    return set(filter(filter_func, objs))


def bibcxx():
    all_compiled_files = (BUILD_DIR / "bibcxx").rglob("*.o")
    bib_objects = set(x.absolute() for x in all_compiled_files)

    bibc_obj_list_path = THIS_DIR / "tmp_bibcxx_objects.txt"

    extra_deps = [
        "esmumps.lib",
        "scotch.lib",
        "scotcherr.lib",
        "metis.lib",
        "medC.lib",
        "medfwrap.lib",
        "hdf5.lib",
        "pthread.lib",
        "TFELSystem.lib",
        "z.lib",
        "MFrontGenericInterface.lib",
        "mkl_intel_lp64_dll.lib",
        "mkl_intel_thread_dll.lib",
        "mkl_core_dll.lib",
        "bibfor.lib",
        "bibfor.exp",
        "bibc.lib",
        "bibc.exp",
    ]
    lib_name = "bibc"
    lib_args = get_default_lib_paths(lib_name)
    args = [
        "LINK.exe",
        f"/IMPLIB:{lib_name}/{lib_name}.lib",
        f"/OUT:{TMP_DIR}/{lib_name}.dll",
        f"@{bibc_obj_list_path}",
    ]

    with open(bibc_obj_list_path, "w") as f:
        f.write("\n".join(map(str, get_default_flags())))
        f.write("\n")
        f.write("\n".join(map(str, lib_args)))
        f.write("\n")
        f.write("\n".join(map(str, bib_objects)))
        f.write("\n")
        f.write("\n".join(map(str, extra_deps)))

    result = call_using_env(args)
    print(result.stdout)
    print(result.stderr)


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
    result = call_using_env(input_args)
    print(result.stdout)
    print(result.stderr)


def bibfor():
    all_bibfor_files = (BUILD_DIR / "bibfor").rglob("*.o")
    bibfor_objects = set(x.absolute() for x in all_bibfor_files)

    cleaned_bibfor_objects = bibfor_objects
    # cleaned_bibfor_objects = filter_objects(bibfor_objects, lambda x: not x.name.startswith("ar_d"))
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
        "pthread.lib",
        "scotch.lib",
        "scotcherr.lib",
        "metis.lib",
        "hdf5.lib",
        "medC.lib",
        "medfwrap.lib",
        "MFrontGenericInterface.lib",
        "mkl_intel_lp64_dll.lib",
        "mkl_intel_thread_dll.lib",
        "mkl_core_dll.lib",
        "bibc.lib",
        "bibcxx.lib",
    ]

    input_args = create_args("bibfor", dep_files, extra_deps, extra_flags=extra_flags)
    input_args += [f"/LIBPATH:{BUILD_DIR}/bibc", f"/LIBPATH:{BUILD_DIR}/bibfor", f"/LIBPATH:{BUILD_DIR}/bibcxx"]

    result = call_using_env(input_args)
    if result.returncode != 0:
        print(" ".join(input_args))
        print(result.stdout)
        print(result.stderr)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Manually Link Code Aster libraries")
    parser.add_argument("--bibc", action="store_true", help="Link the bibc library")
    parser.add_argument("--bibcxx", action="store_true", help="Link the bibcxx library")
    parser.add_argument("--bibfor", action="store_true", help="Link the bibfor library")

    args_ = parser.parse_args()

    if args_.bibc:
        bibc()
    if args_.bibfor:
        bibfor()
    if args_.bibcxx:
        bibcxx()
