import os
import os.path as osp
import pathlib
import shutil

from config import SYMLINK_MAP

SOURCE_DIR = pathlib.Path(os.getenv("CONDA_PREFIX")) / "Library" / "lib/aster"
LIB_DIR = pathlib.Path(os.getenv("CONDA_PREFIX")) / "Library" / "lib"
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


def main(copy_files: bool):
    extlib = ".pyd"
    mods = ["aster", "aster_core", "aster_fonctions", "med_aster", "libaster"]

    py_modules = ["code_aster", "run_aster"]
    for pymod in py_modules:
        ca_module_dir_src = SOURCE_DIR / pymod
        ca_module_dir = SP_DIR / pymod
        if not ca_module_dir.exists():
            shutil.copytree(ca_module_dir_src, ca_module_dir)

    for submodule in mods:
        src = osp.join(DLL_DIR, submodule + extlib)
        dst_dll = SYMLINK_MAP.get(submodule)
        dst = osp.join(SOURCE_DIR, dst_dll)
        if copy_files:
            print(f"Copying {dst} to {src}")
            shutil.copy(dst, src)
        else:
            result = create_symlink(dst, src)
            if not result:
                raise ValueError(f"Failed to create symlink: {src} -> {dst}")


def cli():
    import argparse

    parser = argparse.ArgumentParser(description="Manually Link Code Aster libraries")
    parser.add_argument(
        "--copy",
        action="store_true",
        help="Create symlinks for the Code Aster libraries. Alternatively, copy the files.",
    )
    args = parser.parse_args()
    main(args.copy)


def manual():
    main(False)


if __name__ == "__main__":
    # if no args, use manual mode
    # manual()
    cli()
