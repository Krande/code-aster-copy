import pathlib
import platform
import os.path as osp

import os
import shutil

from conda.config import TMP_DIR

SOURCE_DIR = pathlib.Path(os.getenv("CONDA_PREFIX")) / "Library" / "lib/aster"
#SOURCE_DIR = TMP_DIR
LIB_DIR = pathlib.Path(os.getenv("CONDA_PREFIX")) / "Library" / "lib"
DLL_DIR = pathlib.Path(os.getenv("CONDA_PREFIX")) / "DLLs"
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


def main():
    extlib = ".so"
    libs = ["aster", "bibc", "bibcxx", "bibfor", "AsterMFrOfficialDebug"]
    mods = ["aster", "aster_core", "aster_fonctions", "med_aster"]
    libaster = "libaster.so"
    if platform.system() == "Windows":
        extlib = ".pyd"
        mods.append("libaster")
        libaster = "aster.dll"

    py_modules = ["code_aster", "run_aster"]
    for pymod in py_modules:
        ca_module_dir_src = SOURCE_DIR / pymod
        ca_module_dir = SP_DIR / pymod
        if not ca_module_dir.exists():
            shutil.copytree(ca_module_dir_src, ca_module_dir)

    for lib in libs:
        dll_src = (SOURCE_DIR / lib).with_suffix(".dll")
        lib_src = (SOURCE_DIR / lib).with_suffix(".lib")
        dll_dst = BIN_DIR / dll_src.name
        lib_dst = LIB_DIR / lib_src.name
        print(f"Copying {dll_src} to {dll_dst}")
        shutil.copy(dll_src, dll_dst)
        print(f"Copying {lib_src} to {lib_dst}")
        shutil.copy(lib_src, lib_dst)

    for submodule in mods:
        src = osp.join(DLL_DIR, submodule + extlib)
        dst = osp.join(SOURCE_DIR, libaster)
        # result = create_symlink(dst, src)
        print(f"Copying {dst} to {src}")
        shutil.copy(dst, src)
        # self.symlink_as(src, libaster)


if __name__ == "__main__":
    main()
