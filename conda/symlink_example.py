import pathlib
import platform
import os.path as osp

import os
import shutil

ASTERLIBDIR = pathlib.Path(os.getenv('CONDA_PREFIX')) / "Library" / "lib/aster"
DLL_DIR = pathlib.Path(os.getenv('CONDA_PREFIX')) / "DLLs"


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
    mods = ["aster", "aster_core", "aster_fonctions", "med_aster"]
    libaster = "libaster.so"
    if platform.system() == "Windows":
        extlib = ".pyd"
        mods.append("libaster")
        libaster = "aster.dll"

    for submodule in mods:
        src = osp.join(DLL_DIR, submodule + extlib)
        dst = osp.join(ASTERLIBDIR, libaster)
        # result = create_symlink(dst, src)
        shutil.copy(dst, src)
        # self.symlink_as(src, libaster)


if __name__ == '__main__':
    main()
