# Based on https://learn.microsoft.com/en-us/previous-versions/visualstudio/visual-studio-2008/kkt2hd12(v=vs.90)?redirectedfrom=MSDN
# and https://stackoverflow.com/questions/362830/circular-dependencies-between-dlls-with-visual-studio

import os
import pathlib
import shutil
from enum import Enum

from dotenv import load_dotenv

from msvc_utils import call_using_env

load_dotenv()

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent
THIS_DIR = pathlib.Path(__file__).resolve().parent
BUILD_DIR = ROOT_DIR / "build" / "std" / "debug"
TMP_DIR = THIS_DIR / "temp"
TMP_DIR.mkdir(exist_ok=True)

CONDA_PREFIX_DIR = pathlib.Path(os.getenv("CONDA_PREFIX"))


class CAMod(str, Enum):
    BIBC = "bibc"
    BIBFOR = "bibfor"
    BIBCXX = "bibcxx"
    BIBASTER = "aster"
    ALL = "all"


def run_link(lib_name: CAMod):
    module_dir = BUILD_DIR / str(lib_name.value)
    if lib_name == CAMod.BIBASTER:
        module_dir.mkdir(parents=True, exist_ok=True)
        files = [BUILD_DIR / "bibc" / "supervis" / "python.c.1.o"]
    else:
        files = list(module_dir.rglob("*.o"))

    txt_file = TMP_DIR / f"{lib_name.value}_ofiles.txt"
    out_file = TMP_DIR / f"{lib_name.value}.lib"

    if len(files) == 0:
        raise ValueError(f"No files found in {module_dir}")

    with open(txt_file, "w") as f:
        f.write("\n".join(map(str, files)))

    cmd = [
        "LIB.exe",
        f"/LIBPATH:{CONDA_PREFIX_DIR}/libs",
        f"@{txt_file.as_posix()}",
        f"/OUT:{out_file.as_posix()}",
        "/MACHINE:X64",
        "/VERBOSE"
    ]
    print(" ".join(cmd))

    result = call_using_env(cmd)
    if result.returncode != 0:
        raise ValueError(f"Error linking {lib_name} due to {result.stderr=}, {result.stdout=}")


def main(lib_option: CAMod):
    for lib_name in [CAMod.BIBFOR, CAMod.BIBCXX, CAMod.BIBC, CAMod.BIBASTER]:
        if lib_option != CAMod.ALL and lib_name != lib_option:
            continue
        run_link(lib_name)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Manually Link Code Aster libraries")
    parser.add_argument("--bibc", action="store_true", help="Create .lib for bibc")
    parser.add_argument("--bibfor", action="store_true", help="Create .lib for bibfor")
    parser.add_argument("--bibcxx", action="store_true", help="Create .lib for bibcxx")
    parser.add_argument("--bibaster", action="store_true", help="Create .lib for bibaster")
    parser.add_argument("--all", action="store_true", help="Create .lib for all")

    args = parser.parse_args()

    if not any([args.bibc, args.bibfor, args.bibcxx, args.bibaster, args.all]):
        parser.error("No action requested, add --bibc, --bibfor, --bibcxx, or --all")

    if args.all:
        main(CAMod.ALL)
    elif args.bibc:
        main(CAMod.BIBC)
    elif args.bibfor:
        main(CAMod.BIBFOR)
    elif args.bibcxx:
        main(CAMod.BIBCXX)
    elif args.bibaster:
        main(CAMod.BIBASTER)
    else:
        raise ValueError("Invalid state reached")
