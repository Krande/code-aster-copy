import pathlib
from typing import Iterable

from config import (
    get_obj_list_path,
    get_bibc_compile_files,
    get_bibcxx_compile_files,
    get_bibfor_compile_files,
    CONDA_PREFIX_DIR,
    TMP_DIR,
    ROOT_DIR,
    get_bibaster_compile_files,
)
from msvc_utils import call_using_env


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
        f"/LIBPATH:{CONDA_PREFIX_DIR}/libs",
        f"/LIBPATH:{CONDA_PREFIX_DIR}/include",
        f"/LIBPATH:{ROOT_DIR}/{lib_name}/include",
        f"/LIBPATH:{CONDA_PREFIX_DIR}/Library/lib",
        f"/LIBPATH:{CONDA_PREFIX_DIR}/Library/bin",
        f"/LIBPATH:{TMP_DIR}",
    ]


def run_link(lib_name: str, bib_objects: Iterable[pathlib.Path], extra_deps: list[str]):
    bib_obj_list_path = get_obj_list_path(lib_name)
    lib_args = get_default_lib_paths(lib_name)
    args = [
        "LINK.exe",
        f"/IMPLIB:{TMP_DIR}/{lib_name}.lib",
        f"/OUT:{TMP_DIR}/{lib_name}.dll",
        f"@{bib_obj_list_path}",
    ]
    bib_obj_list_path.parent.mkdir(exist_ok=True, parents=True)
    with open(bib_obj_list_path, "w") as f:
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


def bibaster(use_def: bool = False):
    extra_deps = [
        "esmumps.lib",
        "smumps.lib",
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
        "_raw_bibfor_nodef.lib",
        "_raw_bibc_nodef.lib",
        "_raw_bibcxx_nodef.lib",
    ]
    if use_def:
        extra_deps.append(f"/DEF:{TMP_DIR}/aster_sym.def")

    run_link(lib_name="aster", bib_objects=get_bibaster_compile_files(), extra_deps=extra_deps)


def bibcxx(use_def: bool = False):
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
        "_raw_aster_nodef.lib",
        "_raw_bibfor_nodef.lib",
        "_raw_bibc_nodef.lib",
        "/FORCE:MULTIPLE",
    ]
    if use_def:
        extra_deps.append(f"/DEF:{TMP_DIR}/bibcxx_sym.def")

    run_link(lib_name="bibcxx", bib_objects=get_bibcxx_compile_files(), extra_deps=extra_deps)


def bibc(use_def: bool = False):
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
        "_raw_bibfor_nodef.lib",
        "_raw_bibcxx_nodef.lib",
    ]
    if use_def:
        extra_deps.append(f"/DEF:{TMP_DIR}/bibc_sym.def")

    run_link(lib_name="bibc", bib_objects=get_bibc_compile_files(), extra_deps=extra_deps)


def bibfor(use_def: bool = False):
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
        "_raw_bibc_nodef.lib",
        "_raw_bibcxx_nodef.lib",
    ]
    if use_def:
        extra_deps.append(f"/DEF:{TMP_DIR}/bibfor_sym.def")

    run_link(lib_name="bibfor", bib_objects=get_bibfor_compile_files(), extra_deps=extra_deps)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Manually Link Code Aster libraries")
    parser.add_argument("--bibc", action="store_true", help="Link the bibc library")
    parser.add_argument("--bibcxx", action="store_true", help="Link the bibcxx library")
    parser.add_argument("--bibfor", action="store_true", help="Link the bibfor library")
    parser.add_argument("--bibaster", action="store_true", help="Link the bibaster library")
    parser.add_argument("--use-def", action="store_true", help="Use a .def file to export symbols")

    args_ = parser.parse_args()

    if args_.bibc:
        bibc(args_.use_def)
    if args_.bibfor:
        bibfor(args_.use_def)
    if args_.bibcxx:
        bibcxx(args_.use_def)
    if args_.bibaster:
        bibaster(args_.use_def)
