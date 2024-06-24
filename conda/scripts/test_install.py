# must install dlltracer using p
import os
import pathlib
import shutil

import sys
import ctypes

from config import ROOT_DIR

CONDA_PREFIX = pathlib.Path(os.environ["CONDA_PREFIX"])
ASTEST_DIR = ROOT_DIR / "astest"
BIN_DIR = CONDA_PREFIX / "Library" / "bin"


def init_env():
    ASTER_DIR = CONDA_PREFIX / "Library" / "lib" / "aster"

    os.environ["ASTER_DATADIR"] = (CONDA_PREFIX / "Library" / "share" / "aster").as_posix()
    os.environ["ASTER_LIBDIR"] = ASTER_DIR.as_posix()
    os.environ["ASTER_LOCALEDIR"] = (CONDA_PREFIX / "Library" / "share" / "locale" / "aster").as_posix()
    os.environ["ASTER_ELEMENTSDIR"] = ASTER_DIR.as_posix()


def openmp_debugging():
    import os
    # Set these
    # set KMP_VERSION=1
    # set KMP_AFFINITY=verbose,none
    # set KMP_SETTINGS=1
    # set OMP_DISPLAY_ENV=TRUE
    # set MKL_VERBOSE=1
    # set MKL_DEBUG_CPU_TYPE=5
    # set OMP_NUM_THREADS=1
    # set MKL_NUM_THREADS=1
    # set MKL_DYNAMIC=FALSE
    # set MKL_THREADING_LAYER=INTEL
    os.environ['KMP_VERSION'] = '1'
    os.environ['KMP_AFFINITY'] = 'verbose,none'
    os.environ['KMP_SETTINGS'] = '1'
    os.environ['OMP_DISPLAY_ENV'] = 'TRUE'
    os.environ['MKL_VERBOSE'] = '1'
    os.environ['MKL_DEBUG_CPU_TYPE'] = '5'
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['MKL_DYNAMIC'] = 'FALSE'
    os.environ['MKL_THREADING_LAYER'] = 'INTEL'
    num_threads = os.getenv('OMP_NUM_THREADS')


init_env()

init_str = """
from math import *

import code_aster
from code_aster.Commands import *
from code_aster import CA

"""


# run_ctest --resutest=temp\seq-debug -L submit -L sequential -LE need_data --timefactor=5.0 --only-failed-results -j 4
# run_ctest --resutest=temp/seq -L submit -L sequential -LE need_data --timefactor=5.0 --only-failed-results -j 1
# run_ctest --resutest=temp\mpi -L submit -LE need_data --timefactor=5.0 --only-failed-results -j 6

def run_specific_test(test_name: str, debug_openmp=False):
    comm_file = ASTEST_DIR / f"{test_name}.comm"
    export_file = ASTEST_DIR / f"{test_name}.export"
    mmed_file = ASTEST_DIR / f"{test_name}.mmed"
    mmed_num = 20
    datg_file = ASTEST_DIR / f"{test_name}.datg"
    datg_num = 16

    if export_file.exists():
        export = export_file.read_text(encoding="utf-8")
        for line in export.split("\n"):
            if not line.startswith("F"):
                continue
            line_split = line.split(" ")
            if line_split[1] == "mmed":
                mmed_file = ASTEST_DIR / line_split[2]
                mmed_num = int(line_split[-1])
            if line_split[1] == "datg":
                datg_file = ASTEST_DIR / line_split[2]
                datg_num = int(line_split[-1])

    if mmed_file.exists():
        shutil.copy(mmed_file, f"fort.{mmed_num}")
    if datg_file.exists():
        shutil.copy(datg_file, f"fort.{datg_num}")

    test_file = init_str + comm_file.read_text()

    with open("test_file.py", "w") as f:
        f.write("from test_install import init_env")
        if debug_openmp:
            f.write(", openmp_debugging\n")
            f.write("openmp_debugging()\n")
        else:
            f.write("\n")
        f.write("init_env()\n")
        f.write(test_file)

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
    import dlltracer
    import ctypes
    # ctypes.CDLL(rf"{os.getenv('CONDA_PREFIX')}\Library\lib\aster\bibcxx.dll")
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
    # run_specific_test('adlv100a')
    # run_specific_test('adlv100p')
    # run_specific_test('ahlv100a')
    # run_specific_test('ahlv100t')
    # run_specific_test('hsna106a')  # related to verification of <F> <MODELISA_1> Mesh file issue
    # run_specific_test('zzzz111a')
    run_specific_test('zzzz395s')


if __name__ == "__main__":
    cli()
    # manual()
