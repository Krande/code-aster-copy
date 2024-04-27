import os
import pathlib
import subprocess

from dotenv import load_dotenv

from msvc_utils import call_using_env

load_dotenv()

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent
THIS_DIR = pathlib.Path(__file__).resolve().parent
BUILD_DIR = ROOT_DIR / "build" / "std" / "debug"

CONDA_PREFIX_DIR = pathlib.Path(os.getenv("CONDA_PREFIX"))


def compile_bibfor(specific_file: str = None):
    args = [
        "ifx.exe",
        "/fpp",
        "/MD",
        "/names:lowercase",
        "/assume:underscore",
        "/assume:nobscc",
        "/nologo",
        "/nogen-interfaces",
        "/fpe:0",
        "/traceback",
        "/real-size:64",
        "/double-size:64",
        "/4R8",
        "/4I8",
        "/integer-size:64",
        "/debug:full",
        "/Od",
        "/Qopenmp",
        f"/I{CONDA_PREFIX_DIR}\\Library\\include\\mumps_seq",
        "/I.",
        f"/I{BUILD_DIR}",
        f"/I{CONDA_PREFIX_DIR}\\include",
        f"/I{ROOT_DIR.as_posix()}\\bibfor\\include",
        f"/I{CONDA_PREFIX_DIR}\\Library\\include",
        f"/I{CONDA_PREFIX_DIR}\\Lib\\site-packages\\numpy\\core\\include",
        f"/I{CONDA_PREFIX_DIR}\\libs",
        "/DH5_BUILT_AS_DYNAMIC_LIB",
        "/DASTER_HAVE_INTEL_IFORT",
        "/c",
        "/o",
    ]
    for ffile in (ROOT_DIR / "bibfor").rglob("*.F90"):
        if specific_file and specific_file not in ffile.name:
            continue

        ffile_rel = ffile.relative_to(ROOT_DIR)
        output_ffile = (BUILD_DIR / ffile_rel).with_suffix(".F90.1.o")
        if output_ffile.exists():
            continue
        copy_args = args.copy()
        copy_args.extend([str(output_ffile), str(ffile)])
        print(" ".join(copy_args))
        subprocess.run(copy_args, shell=True, cwd=ROOT_DIR)
        call_using_env(copy_args)
        break


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Manually Link Code Aster libraries")
    parser.add_argument("--file-name", type=str, help="Specific file to compile")

    args = parser.parse_args()
    compile_bibfor(args.file_name)
