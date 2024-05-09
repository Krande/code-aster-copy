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

os.environ["ASTER_DATADIR"] = (CONDA_PREFIX / "Library" / "share" / "aster").as_posix()
os.environ["ASTER_LIBDIR"] = ASTER_DIR.as_posix()
os.environ["ASTER_LOCALEDIR"] = (CONDA_PREFIX / "Library" / "share" / "locale" / "aster").as_posix()
os.environ["ASTER_ELEMENTSDIR"] = ASTER_DIR.as_posix()


def simple():
    from code_aster import CA

    CA.close()


def individual_test():
    pre_reqs = ["medC.dll", "medfwrap.dll"]
    for pre_req in pre_reqs:
        pre_req_path = BIN_DIR / pre_req
        if not pre_req_path.exists():
            raise FileNotFoundError(f"{pre_req_path} does not exist")
        ctypes.CDLL(pre_req_path.as_posix())

    dll_order = ["AsterMFrOfficialDebug.dll", "bibfor.dll", "bibc.dll",  "bibcxx.dll", "aster.dll"]
    for dll_o in dll_order:
        dll = ASTER_DIR / dll_o
        if not dll.exists():
            raise FileNotFoundError(f"{dll} does not exist")
        try:
            ctypes.CDLL(dll.as_posix())
            print(f"Loaded {dll}")
        except OSError as e:
            print(f"Failed to load {dll}: {e}")
            raise e

    import code_aster
    from code_aster import CA


def trace_test():
    import ctypes
    ctypes.CDLL(rf"{os.getenv('CONDA_PREFIX')}\Library\lib\aster\bibcxx.dll")
    sys.path.insert(0, ASTER_DIR.as_posix())
    with dlltracer.Trace(out=sys.stdout):
        import code_aster
        from code_aster import CA


def cli():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--simple", action="store_true")
    parser.add_argument("--trace", action="store_true")
    parser.add_argument("--seq", action="store_true")

    args = parser.parse_args()

    if args.simple:
        simple()

    if args.trace:
        trace_test()

    if args.seq:
        individual_test()


if __name__ == "__main__":
    cli()
