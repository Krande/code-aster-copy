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
    MODS, SYMLINK_MAP,
)
from msvc_utils import call_using_env


# https://learn.microsoft.com/en-us/cpp/build/reference/linker-options?view=msvc-170
def get_default_flags():
    return ["/nologo", "/subsystem:console", "/DLL", "/MACHINE:X64", "/DEBUG"]


def get_default_lib_paths(lib_name: str):
    return [
        f"/LIBPATH:{TMP_DIR}",
        f"/LIBPATH:{CONDA_PREFIX_DIR}/libs",
        f"/LIBPATH:{CONDA_PREFIX_DIR}/include",
        f"/LIBPATH:{ROOT_DIR}/{lib_name}/include",
        f"/LIBPATH:{CONDA_PREFIX_DIR}/Library/lib",
        f"/LIBPATH:{CONDA_PREFIX_DIR}/Library/bin",
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
        f"/IMPLIB:{TMP_DIR}/{lib_name}.lib",
        f"/OUT:{TMP_DIR}/{lib_name}.dll",
    ]

    bib_obj_list_path.parent.mkdir(exist_ok=True, parents=True)
    with open(bib_obj_list_path, "w") as f:
        f.write("\n".join(map(str, args)))
        f.write("\n")
        f.write("\n".join(map(str, get_default_flags())))
        f.write("\n")
        f.write("\n".join(map(str, lib_args)))
        f.write("\n")
        f.write("\n".join(map(str, bib_objects)))
        f.write("\n")
        f.write("\n".join(map(str, extra_deps)))

    result = call_using_env(["LINK.exe", f"@{bib_obj_list_path}"])
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise ValueError(f"Error linking {lib_name}")

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

    extra_deps = core_lib_deps
    extra_deps += SHARED_DEPS + ["/WHOLEARCHIVE:aster.lib"] + ["aster.exp"]
    run_link(lib_name="aster", bib_objects=get_bibaster_compile_files(), extra_deps=extra_deps)


def bibcxx():
    deps = LIB_DEPENDENCIES.get(CAMod.BIBCXX)
    core_lib_deps = eval_deps(deps)

    extra_deps = core_lib_deps
    extra_deps += SHARED_DEPS + ["/WHOLEARCHIVE:bibcxx.lib"] + ["bibcxx.exp"]

    run_link(lib_name="bibcxx", bib_objects=get_bibcxx_compile_files(), extra_deps=extra_deps)


def bibc():
    deps = LIB_DEPENDENCIES.get(CAMod.BIBC)
    core_lib_deps = eval_deps(deps)

    extra_deps = core_lib_deps
    extra_deps += SHARED_DEPS + ["/WHOLEARCHIVE:bibc.lib"] + ["bibc.exp"]

    run_link(lib_name="bibc", bib_objects=get_bibc_compile_files(), extra_deps=extra_deps)


def bibfor():
    deps = LIB_DEPENDENCIES.get(CAMod.BIBFOR)
    core_lib_deps = eval_deps(deps)

    extra_deps = core_lib_deps
    extra_deps += SHARED_DEPS + ["/WHOLEARCHIVE:bibfor.lib"] + ["bibfor.exp"]
    run_link(lib_name="bibfor", bib_objects=get_bibfor_compile_files(), extra_deps=extra_deps)


def aster_lib_complete():
    """Link all the libraries into a single dll"""
    ca_module_deps = LIB_DEPENDENCIES.get(CAMod.ALL)
    lib_deps = eval_deps(ca_module_deps)

    extra_deps = SHARED_DEPS + lib_deps
    # extra_deps += ["aster.exp", "bibc.exp", "bibcxx.exp", "bibfor.exp"]
    extra_deps += [
        # "/WHOLEARCHIVE:aster.lib",
        # "/WHOLEARCHIVE:bibcxx.lib",
        # "/WHOLEARCHIVE:bibc.lib",
        # "/WHOLEARCHIVE:bibfor.lib",
        "/FORCE:MULTIPLE",
    ]

    object_files = get_bibfor_compile_files()
    object_files += get_bibc_compile_files()
    object_files += get_bibcxx_compile_files()
    object_files += get_bibaster_compile_files()

    run_link(lib_name="astermain", bib_objects=object_files, extra_deps=extra_deps)


def cli():
    import argparse

    parser = argparse.ArgumentParser(description="Manually Link Code Aster libraries")
    parser.add_argument("--bibc", action="store_true", help="Link the bibc library")
    parser.add_argument("--bibcxx", action="store_true", help="Link the bibcxx library")
    parser.add_argument("--bibfor", action="store_true", help="Link the bibfor library")
    parser.add_argument("--bibaster", action="store_true", help="Link the bibaster library")
    parser.add_argument("--complete", action="store_true", help="Link all library into a single dll")
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
    if args_.complete:
        aster_lib_complete()


def manual():
    shutil.rmtree(TMP_DIR, ignore_errors=True)
    for mod in [CAMod.BIBCXX, CAMod.BIBFOR, CAMod.BIBC, CAMod.LIBASTER]:
        run_lib(mod, DEFOption.USE_DEF)

    aster_lib_complete()

    # bibfor()
    # bibc()
    # bibcxx()
    # bibaster()

    # Run linking again
    # for mod in MODS:
    #     dll_name = SYMLINK_MAP.get(mod)
    #     shutil.copy(TMP_DIR / dll_name, (TMP_DIR / mod).with_suffix(".pyd"))


if __name__ == "__main__":
    # cli()
    manual()
