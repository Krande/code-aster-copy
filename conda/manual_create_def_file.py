import pathlib

from config import TMP_DIR, LIB_RAW_PREFIX
from manual_lib import run_lib
from msvc_utils import call_using_env


def create_symbols_file(module_name: str) -> None:
    symbol_output_file = (TMP_DIR / f"{LIB_RAW_PREFIX}{module_name}_sym").with_suffix(".txt")
    bib_lib_file = symbol_output_file.with_name(f"{LIB_RAW_PREFIX}{module_name}.lib")

    if not bib_lib_file.exists():
        bib_lib_file = run_lib(module_name, use_def=False)

    result = call_using_env(["dumpbin", "/symbols", bib_lib_file, ">", symbol_output_file])
    if result.returncode != 0:
        raise ValueError(f"Error dumping symbols for {module_name} due to {result.stderr=}, {result.stdout=}")


def symbol_out_to_def_file(symbol_file: pathlib.Path):
    stem_clean = symbol_file.stem.replace(LIB_RAW_PREFIX, "")
    def_file = symbol_file.with_name(f"{stem_clean}.def")
    written_symbols = set()
    with open(def_file, "w") as fdef:
        fdef.write("EXPORTS\n")
        with open(symbol_file) as f:
            for line in f:
                if "UNDEF" in line:
                    continue
                if "SECT1" not in line:
                    continue
                if "()" not in line:
                    continue
                if "00000000" not in line:
                    continue
                symbol_name = line.split("|")[-1].strip()
                if symbol_name in written_symbols:
                    continue
                written_symbols.add(symbol_name)
                fdef.write(f"{symbol_name}\n")


def create_symbol_def(module_name: str):
    symbol_output_file = (TMP_DIR / f"{LIB_RAW_PREFIX}{module_name}_sym").with_suffix(".txt")
    if not symbol_output_file.exists():
        create_symbols_file(module_name)

    symbol_out_to_def_file(symbol_output_file)


if __name__ == "__main__":
    for mod in ["bibc", "bibfor", "bibcxx"]:
        create_symbol_def(mod)
