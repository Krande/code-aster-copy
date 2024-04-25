# Based on https://learn.microsoft.com/en-us/previous-versions/visualstudio/visual-studio-2008/kkt2hd12(v=vs.90)?redirectedfrom=MSDN
# and https://stackoverflow.com/questions/362830/circular-dependencies-between-dlls-with-visual-studio

import os
import pathlib
import subprocess
from enum import Enum

from dotenv import load_dotenv

load_dotenv()

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent
THIS_DIR = pathlib.Path(__file__).resolve().parent
BUILD_DIR = ROOT_DIR / "build" / "std" / "debug"

CONDA_PREFIX_DIR = pathlib.Path(os.getenv("CONDA_PREFIX"))


class CAModule(Enum, str):
    BIBC = "bibc"
    BIBFOR = "bibfor"
    BIBCXX = "bibcxx"
    BIBCXXASTER = "bibcxxaster"
    ALL = "all"


def main(ca_module: CAModule):
    for lib_name in [CAModule.BIBFOR, CAModule.BIBCXX, CAModule.BIBC, CAModule.BIBCXXASTER]:
        module_dir = BUILD_DIR / str(lib_name)
        if lib_name == CAModule.BIBCXXASTER:
            files = [BUILD_DIR / "bibc" / "supervis" / "python.c.2.o"]
        else:
            files = list(module_dir.rglob("*.o"))
        txt_file = module_dir / f"{lib_name}_ofiles.txt"
        out_file = module_dir / f"{lib_name}.lib"

        if ca_module != CAModule.ALL and lib_name != ca_module:
            continue

        if len(files) == 0:
            raise ValueError(f"No files found in {module_dir}")

        with open(txt_file, "w") as f:
            f.write("\n".join(map(str, files)))

        args = [
            "LIB.exe",
            f"/LIBPATH:{CONDA_PREFIX_DIR}/libs",
            f"@{txt_file.as_posix()}",
            f"/OUT:{out_file.as_posix()}",
            "/MACHINE:X64",
        ]
        print(" ".join(args))

        subprocess.run(args, shell=True, cwd=module_dir)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Manually Link Code Aster libraries")
    parser.add_argument("--bibc", action="store_true", help="Create .lib for bibc")
    parser.add_argument("--bibfor", action="store_true", help="Create .lib for bibfor")
    parser.add_argument("--bibcxx", action="store_true", help="Create .lib for bibcxx")
    parser.add_argument("--bibcxxaster", action="store_true", help="Create .lib for bibcxxaster")
    parser.add_argument("--all", action="store_true", help="Create .lib for all")

    args = parser.parse_args()

    if not any([args.bibc, args.bibfor, args.bibcxx, args.bibcxxaster, args.all]):
        parser.error("No action requested, add --bibc, --bibfor, --bibcxx, or --all")

    if args.all:
        main(CAModule.ALL)
    elif args.bibc:
        main(CAModule.BIBC)
    elif args.bibfor:
        main(CAModule.BIBFOR)
    elif args.bibcxx:
        main(CAModule.BIBCXX)
    elif args.bibcxxaster:
        main(CAModule.BIBCXXASTER)
    else:
        raise ValueError("Invalid state reached")
