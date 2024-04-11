# Based on https://learn.microsoft.com/en-us/previous-versions/visualstudio/visual-studio-2008/kkt2hd12(v=vs.90)?redirectedfrom=MSDN
# and https://stackoverflow.com/questions/362830/circular-dependencies-between-dlls-with-visual-studio

import os
import pathlib
import subprocess

from dotenv import load_dotenv

load_dotenv()

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent
THIS_DIR = pathlib.Path(__file__).resolve().parent
BUILD_DIR = ROOT_DIR / "build" / "std" / "debug"


CONDA_PREFIX_DIR = pathlib.Path(os.getenv("CONDA_PREFIX"))


def main(compile_bibc=False, compile_bibfor=False, compile_bibcxx=False):
    for lib_name in ["bibfor", "bibcxx", "bibc"]:
        module_dir = BUILD_DIR / lib_name
        files = list(module_dir.rglob("*.o"))
        txt_file = module_dir / f"{lib_name}_ofiles.txt"
        out_file = module_dir / f"{lib_name}.lib"

        if not compile_bibc and lib_name == "bibc":
            continue

        if not compile_bibfor and lib_name == "bibfor":
            continue

        if not compile_bibcxx and lib_name == "bibcxx":
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
    parser.add_argument("--bibc", action="store_true", help="Create .lib and .exp for the bibc library")
    parser.add_argument("--bibfor", action="store_true", help="Create .lib and .exp for the bibfor library")
    parser.add_argument("--bibcxx", action="store_true", help="Create .lib and .exp for the bibcxx library")

    args = parser.parse_args()

    main(compile_bibc=args.bibc, compile_bibfor=args.bibfor, compile_bibcxx=args.bibcxx)