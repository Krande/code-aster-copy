import os
import os.path as osp
import pathlib
import shutil

from config import TMP_DIR

SOURCE_DIR = pathlib.Path(os.getenv("CONDA_PREFIX")) / "Library" / "lib/aster"
# SOURCE_DIR = TMP_DIR
LIB_DIR = pathlib.Path(os.getenv("CONDA_PREFIX")) / "Library" / "lib"
# DLL_DIR = pathlib.Path(os.getenv("CONDA_PREFIX")) / "DLLs"
DLL_DIR = SOURCE_DIR
BIN_DIR = pathlib.Path(os.getenv("CONDA_PREFIX")) / "Library" / "bin"
SP_DIR = pathlib.Path(os.getenv("CONDA_PREFIX")) / "Lib" / "site-packages"


def create_symlink(source, link_name):
    """
    Attempts to create a symbolic link on Windows.

    Args:
    source (str): The path to the target file or directory.
    link_name (str): The path where the symbolic link should be created.

    Returns:
    bool: True if the symlink was created successfully, False otherwise.
    """
    try:
        os.symlink(source, link_name)
        print(f"Symlink created successfully: {link_name} -> {source}")
        return True
    except OSError as e:
        print(f"Failed to create symlink: {e}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False


def main(use_symlink: bool):
    extlib = ".pyd"
    libs = ["aster", "bibc", "bibcxx", "bibfor", "AsterMFrOfficialDebug"]
    mods = ["aster", "aster_core", "aster_fonctions", "med_aster", "libaster"]
    libaster = "aster.dll"

    py_modules = ["code_aster", "run_aster"]
    for pymod in py_modules:
        ca_module_dir_src = SOURCE_DIR / pymod
        ca_module_dir = SP_DIR / pymod
        if not ca_module_dir.exists():
            shutil.copytree(ca_module_dir_src, ca_module_dir)

    # for lib in libs:
    #     dll_src = (SOURCE_DIR / lib).with_suffix(".dll")
    #     lib_src = (SOURCE_DIR / lib).with_suffix(".lib")
    #     dll_dst = BIN_DIR / dll_src.name
    #     lib_dst = LIB_DIR / lib_src.name
    #     print(f"Copying {dll_src} to {dll_dst}")
    #     shutil.copy(dll_src, dll_dst)
    #     print(f"Copying {lib_src} to {lib_dst}")
    #     shutil.copy(lib_src, lib_dst)

    for submodule in mods:
        src = osp.join(DLL_DIR, submodule + extlib)
        dst = osp.join(SOURCE_DIR, libaster)
        if use_symlink:
            result = create_symlink(dst, src)
            if not result:
                raise ValueError(f"Failed to create symlink: {src} -> {dst}")
        else:
            print(f"Copying {dst} to {src}")
            shutil.copy(dst, src)


def cli():
    import argparse

    parser = argparse.ArgumentParser(description="Manually Link Code Aster libraries")
    parser.add_argument(
        "--symlink",
        action="store_true",
        help="Create symlinks for the Code Aster libraries. Alternatively, copy the files.",
    )
    args = parser.parse_args()
    main(args.symlink)


def manual():
    main(False)


if __name__ == "__main__":
    # if no args, use manual mode
    # manual()
    cli()
