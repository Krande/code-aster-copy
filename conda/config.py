import os
import pathlib
from enum import Enum
from typing import Iterable

from dotenv import load_dotenv

load_dotenv()

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent
THIS_DIR = pathlib.Path(__file__).resolve().parent
BUILD_DIR = ROOT_DIR / "build" / "std" / "debug"
TMP_DIR = THIS_DIR / "temp"
LIB_RAW_PREFIX = ""

CONDA_PREFIX_DIR = pathlib.Path(os.getenv("CONDA_PREFIX"))
MODS = ["aster", "aster_core", "aster_fonctions", "med_aster", "libaster"]


class CAMod(str, Enum):
    BIBC = "bibc"
    BIBFOR = "bibfor"
    BIBCXX = "bibcxx"
    LIBASTER = "aster"
    MFRONT = "mfront"
    ALL = "all"


class DEFOption(str, Enum):
    USE_DEF = "use_def"
    USE_BLANK_DEF = "use_blank_def"
    NO_DEF = "no_def"


SYMLINK_MAP = {
    "libaster": "bibcxx.dll",
    "aster": "bibc.dll",
    "aster_core": "bibc.dll",
    "aster_fonctions": "bibc.dll",
    "med_aster": "bibc.dll",
}

LIB_DEPENDENCIES = {
    CAMod.BIBC: [CAMod.BIBCXX, CAMod.BIBFOR ],
    CAMod.BIBFOR: [CAMod.BIBC, CAMod.BIBCXX ],
    CAMod.BIBCXX: [CAMod.LIBASTER, CAMod.BIBC, CAMod.BIBFOR],
    CAMod.LIBASTER: [CAMod.BIBCXX, CAMod.BIBC, CAMod.BIBFOR],
    CAMod.MFRONT: [CAMod.BIBCXX],
}
# link passes
# 1. make .lib files for each module using /DEF:{def_file} option
# 2. Create bibc.dll from bibc.exp and bibc.lib

DEF_FILE_MAP = {
    CAMod.BIBC: ROOT_DIR / "bibc/bibc.def",
    CAMod.BIBCXX: ROOT_DIR / "bibcxx/bibcxx.def",
    CAMod.BIBFOR: ROOT_DIR / "bibfor/bibfor.def",
    CAMod.LIBASTER: ROOT_DIR / "bibc/aster.def",
}


class CompileStage(str, Enum):
    LIB = "lib"
    LINK = "link"


def get_obj_list_path(lib_name: str, cstage: CompileStage) -> pathlib.Path:
    return TMP_DIR / "inputs" / f"{LIB_RAW_PREFIX}{lib_name}_{cstage.value}.txt"


def get_obj_sym_file(lib_name: str) -> pathlib.Path:
    return (TMP_DIR / f"{LIB_RAW_PREFIX}{lib_name}_sym").with_suffix(".txt")


def get_lib_file(lib_name: str, def_opt: DEFOption) -> pathlib.Path:
    if def_opt == DEFOption.USE_BLANK_DEF:
        lib_file_name = f"{LIB_RAW_PREFIX}{lib_name}"
    elif def_opt == DEFOption.USE_DEF:
        lib_file_name = f"{LIB_RAW_PREFIX}{lib_name}"
    else:  # NO_DEF
        lib_file_name = f"{LIB_RAW_PREFIX}{lib_name}"

    return (TMP_DIR / lib_file_name).with_suffix(".lib")


def get_bibc_compile_files():
    files = list((BUILD_DIR / "bibc").rglob("*.o"))
    # remove the python.o file
    files_clean = set([x for x in files if "python.c" not in x.stem])
    diff_files = set(files).difference(files_clean)
    if len(diff_files) == 0:
        raise ValueError(f"Removed too many files {diff_files}")
    files = files_clean
    return files


def get_bibaster_compile_files() -> Iterable[pathlib.Path]:
    return [BUILD_DIR / "bibc" / "supervis" / "python.c.2.o"]


def get_bibcxx_compile_files():
    files = list((BUILD_DIR / "bibcxx").rglob("*.o"))
    files_clean = set([x for x in files if "MFrontBehaviour.cxx" not in x.stem])
    diff_files = set(files).difference(files_clean)
    if len(diff_files) == 0:
        raise ValueError(f"Removed too many files {diff_files}")
    files = files_clean
    return files


def get_mfront_compile_files():
    return list((BUILD_DIR / "bibcxx" / "MFront").rglob("*.o"))


def get_bibfor_compile_files():
    return list((BUILD_DIR / "bibfor").rglob("*.o"))
