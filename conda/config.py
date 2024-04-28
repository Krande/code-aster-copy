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
LIB_RAW_PREFIX = "_raw_"
TMP_DIR.mkdir(exist_ok=True)

CONDA_PREFIX_DIR = pathlib.Path(os.getenv("CONDA_PREFIX"))


class CAMod(str, Enum):
    BIBC = "bibc"
    BIBFOR = "bibfor"
    BIBCXX = "bibcxx"
    BIBASTER = "aster"
    ALL = "all"


class DEFOption(str, Enum):
    USE_DEF = "use_def"
    USE_BLANK_DEF = "use_blank_def"
    NO_DEF = "no_def"


def get_obj_list_path(lib_name: str) -> pathlib.Path:
    return TMP_DIR / "inputs" / f"{LIB_RAW_PREFIX}{lib_name}_ofiles.txt"


def get_obj_sym_file(lib_name: str) -> pathlib.Path:
    return (TMP_DIR / f"{LIB_RAW_PREFIX}{lib_name}_sym").with_suffix(".txt")


def get_lib_file(lib_name: str, def_opt: DEFOption) -> pathlib.Path:
    if def_opt == DEFOption.USE_BLANK_DEF:
        lib_file_name = f"{LIB_RAW_PREFIX}{lib_name}_blankdef"
    elif def_opt == DEFOption.USE_DEF:
        lib_file_name = f"{LIB_RAW_PREFIX}{lib_name}_def"
    else:  # NO_DEF
        lib_file_name = f"{LIB_RAW_PREFIX}{lib_name}_nodef"

    return (TMP_DIR / lib_file_name).with_suffix(".lib")


def get_bibc_compile_files(skip_pythonc=True):
    files = list((BUILD_DIR / "bibc").rglob("*.o"))
    if skip_pythonc:
        # remove the python.o file
        files_clean = [x for x in files if "python.c.2.o" not in str(x)]
        if len(files) - len(files_clean) != 1:
            raise ValueError(f"python.c.2.o not found in {files}")
        files = files_clean
    return files


def get_bibaster_compile_files() -> Iterable[pathlib.Path]:
    return [BUILD_DIR / "bibc" / "supervis" / "python.c.2.o"]


def get_bibcxx_compile_files(skip_certain_files=False):
    files = list((BUILD_DIR / "bibcxx").rglob("*.o"))
    if skip_certain_files:
        to_be_removed = {
            "ConstantFieldOnCells.cxx.2.o",
            "ElementaryTerm.cxx.2.o",
            "FieldOnCells.cxx.2.o",
            "BehaviourDefinition.cxx.2.o",
            "ElementaryModeling.cxx.2.o",
            "FieldOnNodes.cxx.2.o",
        }
        files_cleaned = set([x for x in files if not any(y in str(x) for y in to_be_removed)])
        result = files_cleaned.intersection(to_be_removed)
        if len(result) > 0:
            raise ValueError(f"These files where not removed {result} from source")
        files = files_cleaned

    return files


def get_bibfor_compile_files():
    return (BUILD_DIR / "bibfor").rglob("*.o")
