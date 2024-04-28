import os
import pathlib
from enum import Enum

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
