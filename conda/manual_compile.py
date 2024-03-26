import os
import subprocess
import pathlib
from dotenv import load_dotenv

load_dotenv()

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent
THIS_DIR = pathlib.Path(__file__).resolve().parent
BUILD_DIR = ROOT_DIR / "build" / "debug"

CONDA_PREFIX_DIR = pathlib.Path(os.getenv("CONDA_PREFIX"))


def compile_bibfor():
    args = [
        (ROOT_DIR / "call_compile.bat").as_posix(),
        "ifx.exe",
        "-fpp",
        "/MD",
        "/names:lowercase",
        "/nologo",
        "/nogen-interfaces",
        "/fpe:0",
        "/traceback",
        "/real-size:64",
        "/double-size:64",
        "/4R8",
        "/4I8",
        "/integer-size:64",
        "/QxCORE-AVX2",
        "/debug:full",
        "/Od",
        "/Qopenmp",
        f"/module:{BUILD_DIR.as_posix()}",
        f"/I{CONDA_PREFIX_DIR}\\Library\\include\\mumps_seq",
        "/I.",
        f"/I{BUILD_DIR}",
        f"/I{CONDA_PREFIX_DIR}\\include",
        f"/I{ROOT_DIR.as_posix()}/bibfor/include",
        f"/I{CONDA_PREFIX_DIR}\\Library\\include",
        f"/I:{CONDA_PREFIX_DIR}\\Lib\\site-packages\\numpy\\core\\include",
        f"/I{CONDA_PREFIX_DIR}\\libs",
        "/DH5_BUILT_AS_DYNAMIC_LIB",
        "/DASTER_HAVE_INTEL_IFORT",
        "/c",
        "/o",
    ]
    for ffile in (ROOT_DIR / "bibfor").rglob("*.F90"):
        ffile_rel = ffile.relative_to(ROOT_DIR)
        output_ffile = (BUILD_DIR / ffile_rel).with_suffix(".F90.1.o")
        if output_ffile.exists():
            continue
        copy_args = args.copy()
        copy_args.extend([output_ffile.as_posix(), ffile.as_posix()])
        print(" ".join(copy_args))
        subprocess.run(copy_args, shell=True, cwd=BUILD_DIR)
        break


def main():
    compile_bibfor()


if __name__ == "__main__":
    main()
