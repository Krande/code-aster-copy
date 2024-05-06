# must install dlltracer using p
import os
import pathlib

import dlltracer
import sys
import ctypes

from config import TMP_DIR

CONDA_PREFIX = pathlib.Path(os.environ["CONDA_PREFIX"])
ASTER_DIR = CONDA_PREFIX / "Library" / "lib" / "aster"
BIN_DIR = CONDA_PREFIX / "Library" / "bin"


def main():
    pre_reqs = ["medC.dll", "medfwrap.dll"]
    for pre_req in pre_reqs:
        pre_req_path = BIN_DIR / pre_req
        if not pre_req_path.exists():
            raise FileNotFoundError(f"{pre_req_path} does not exist")
        ctypes.CDLL(pre_req_path.as_posix())

    dll_order = ["AsterMFrOfficialDebug.dll", "aster.dll", "bibc.dll", "bibfor.dll", "bibcxx.dll"]
    for dll_o in dll_order:
        dll = ASTER_DIR / dll_o
        if not dll.exists():
            raise FileNotFoundError(f"{dll} does not exist")
        try:
            ctypes.CDLL(dll.as_posix())
        except OSError as e:
            print(f"Failed to load {dll}: {e}")
            raise e

    # sys.path.insert(0, TMP_DIR.as_posix())
    sys.path.insert(0, ASTER_DIR.as_posix())
    with dlltracer.Trace(out=sys.stdout):
        import code_aster
        from code_aster import CA


if __name__ == "__main__":
    main()
