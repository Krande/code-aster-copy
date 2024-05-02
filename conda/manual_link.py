import pathlib
import shutil
from typing import Iterable

from conda.manual_lib import run_lib
from config import (
    get_obj_list_path,
    get_bibc_compile_files,
    get_bibcxx_compile_files,
    get_bibfor_compile_files,
    CONDA_PREFIX_DIR,
    TMP_DIR,
    ROOT_DIR,
    get_bibaster_compile_files,
    LIB_RAW_PREFIX,
    CAMod,
    DEFOption,
    LIB_DEPENDENCIES,
    CompileStage,
    MODS,
)
from msvc_utils import call_using_env


# https://learn.microsoft.com/en-us/cpp/build/reference/linker-options?view=msvc-170
def get_default_flags():
    return ["/nologo", "/MANIFEST", "/subsystem:console", "/DLL"]


def get_default_lib_paths(lib_name: str):
    return [
        f"/LIBPATH:{CONDA_PREFIX_DIR}/libs",
        f"/LIBPATH:{CONDA_PREFIX_DIR}/include",
        f"/LIBPATH:{ROOT_DIR}/{lib_name}/include",
        f"/LIBPATH:{CONDA_PREFIX_DIR}/Library/lib",
        f"/LIBPATH:{CONDA_PREFIX_DIR}/Library/bin",
        f"/LIBPATH:{TMP_DIR}",
    ]


SHARED_DEPS = [
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
]


def run_link(lib_name: str, bib_objects: Iterable[pathlib.Path], extra_deps: list[str]):
    bib_obj_list_path = get_obj_list_path(lib_name, CompileStage.LINK)
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
    return result


def eval_deps(deps):
    """The pre-check if libs exists or not can be done here."""
    prefix = LIB_RAW_PREFIX
    core_lib_deps = []
    for i, dep_enum in enumerate(deps):
        dep = f"{prefix}{dep_enum.value}.lib"
        lib_path = TMP_DIR / dep
        if not lib_path.exists():
            run_lib(dep_enum, def_opt=DEFOption.USE_DEF)
        dep_core_path = (TMP_DIR / dep_enum.value).with_suffix(".lib")
        exp_file = lib_path.with_suffix(".exp")
        # if dep_core_path.with_suffix(".exp").exists():
        #     # This means that /DEF has been applied on this library
        #     core_lib_deps.extend([dep_core_path.name])
        #     continue
        # else:
        core_lib_deps.extend([dep])


    return core_lib_deps


def bibaster():
    deps = LIB_DEPENDENCIES.get(CAMod.LIBASTER)
    core_lib_deps = eval_deps(deps)

    extra_deps = SHARED_DEPS + core_lib_deps

    result = run_link(lib_name="aster", bib_objects=get_bibaster_compile_files(), extra_deps=extra_deps)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise ValueError("Error linking aster")


def bibcxx():
    deps = LIB_DEPENDENCIES.get(CAMod.BIBCXX)
    core_lib_deps = eval_deps(deps)

    extra_deps = SHARED_DEPS + core_lib_deps

    result = run_link(lib_name="bibcxx", bib_objects=get_bibcxx_compile_files(), extra_deps=extra_deps)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise ValueError("Error linking bibcxx")


def bibc():
    deps = LIB_DEPENDENCIES.get(CAMod.BIBC)
    core_lib_deps = eval_deps(deps)
    extra_deps = SHARED_DEPS + core_lib_deps

    result = run_link(lib_name="bibc", bib_objects=get_bibc_compile_files(), extra_deps=extra_deps)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise ValueError("Error linking bibc")


def bibfor():
    deps = LIB_DEPENDENCIES.get(CAMod.BIBFOR)
    core_lib_deps = eval_deps(deps)

    extra_deps = SHARED_DEPS + core_lib_deps

    result = run_link(lib_name="bibfor", bib_objects=get_bibfor_compile_files(), extra_deps=extra_deps)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise ValueError("Error linking bibfor")


def cli():
    import argparse

    parser = argparse.ArgumentParser(description="Manually Link Code Aster libraries")
    parser.add_argument("--bibc", action="store_true", help="Link the bibc library")
    parser.add_argument("--bibcxx", action="store_true", help="Link the bibcxx library")
    parser.add_argument("--bibfor", action="store_true", help="Link the bibfor library")
    parser.add_argument("--bibaster", action="store_true", help="Link the bibaster library")
    parser.add_argument("--use-def", action="store_true", help="Use a .def file to export symbols")

    args_ = parser.parse_args()

    if args_.bibc:
        bibc()
    if args_.bibfor:
        bibfor()
    if args_.bibcxx:
        bibcxx()
    if args_.bibaster:
        bibaster()


def manual():
    shutil.rmtree(TMP_DIR, ignore_errors=True)
    bibc()
    bibfor()
    bibcxx()
    bibaster()
    # Run linking again
    for mod in MODS:
        shutil.copy(TMP_DIR / "aster.dll", (TMP_DIR / mod).with_suffix(".pyd"))


if __name__ == "__main__":
    # cli()
    manual()
