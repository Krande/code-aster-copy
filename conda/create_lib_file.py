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
BIBFOR_BUILD_DIR = BUILD_DIR / "bibfor"

CONDA_PREFIX_DIR = pathlib.Path(os.getenv("CONDA_PREFIX"))


def main():
    files = list(BIBFOR_BUILD_DIR.rglob("*.F90.1.o"))
    txt_file = BIBFOR_BUILD_DIR / "bibfor_ofiles.txt"
    out_file = BIBFOR_BUILD_DIR / "bibfor.lib"

    if len(files) == 0:
        raise ValueError(f"No files found in {BIBFOR_BUILD_DIR}")

    with open(txt_file, "w") as f:
        f.write("\n".join(map(str, files)))

    # "lib /DEF @conda/bibfor_ofiles.txt /OUT:bibfor.lib"
    args = [
        "lib",
        "/DEF",
        f"@{txt_file.as_posix()}",
        f"/OUT:{out_file.as_posix()}",
    ]
    print(" ".join(args))

    subprocess.run(args, shell=True, cwd=BIBFOR_BUILD_DIR)


if __name__ == "__main__":
    main()
