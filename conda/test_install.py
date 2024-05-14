# must install dlltracer using p
import os
import pathlib
import shutil

import dlltracer
import sys
import ctypes

from config import ROOT_DIR

CONDA_PREFIX = pathlib.Path(os.environ["CONDA_PREFIX"])
ASTER_DIR = CONDA_PREFIX / "Library" / "lib" / "aster"
ASTEST_DIR = ROOT_DIR / "astest"
BIN_DIR = CONDA_PREFIX / "Library" / "bin"

os.environ["ASTER_DATADIR"] = (CONDA_PREFIX / "Library" / "share" / "aster").as_posix()
os.environ["ASTER_LIBDIR"] = ASTER_DIR.as_posix()
os.environ["ASTER_LOCALEDIR"] = (CONDA_PREFIX / "Library" / "share" / "locale" / "aster").as_posix()
os.environ["ASTER_ELEMENTSDIR"] = ASTER_DIR.as_posix()

init_str = """
from math import *

import code_aster
from code_aster.Commands import *
from code_aster import CA

"""


def run_specific_test(test_name: str):
    comm_file = ASTEST_DIR / f"{test_name}.comm"
    mmed_file = ASTEST_DIR / f"{test_name}.mmed"
    if mmed_file.exists():
        shutil.copy(mmed_file, "fort.20")

    test_file = init_str + comm_file.read_text()
    exec(test_file)


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

    dll_order = ["AsterMFrOfficialDebug.dll", "bibfor.dll", "bibc.dll", "bibcxx.dll", "aster.dll"]
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
    parser.add_argument("--test-file", type=str, help="Run a specific test file")
    parser.add_argument("--manual", action="store_true")

    args = parser.parse_args()

    if args.simple:
        simple()

    if args.trace:
        trace_test()

    if args.seq:
        individual_test()

    if args.test_file:
        run_specific_test(args.test_file)

    if args.manual:
        manual()


def manual():
    # run_specific_test("comp010i")
    run_specific_test('adlv100a')


if __name__ == "__main__":
    cli()
    # manual()
