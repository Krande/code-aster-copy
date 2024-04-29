import pathlib

from config import TMP_DIR, LIB_RAW_PREFIX, get_obj_sym_file, get_lib_file, DEFOption
from manual_lib import run_lib
from msvc_utils import call_using_env


def create_symbols_file(module_name: str, def_option: DEFOption) -> None:
    symbol_output_file = get_obj_sym_file(module_name)
    bib_lib_file = get_lib_file(module_name, def_option)

    if not bib_lib_file.exists():
        run_lib(module_name, def_opt=def_option)

    result = call_using_env(["dumpbin", "/symbols", bib_lib_file, ">", symbol_output_file])
    if result.returncode != 0:
        raise ValueError(f"Error dumping symbols for {module_name} due to {result.stderr=}, {result.stdout=}")


def iter_symbol_names(symbol_file: pathlib.Path | str):
    written_symbols = set()
    with open(symbol_file) as f:
        for line in f:
            if "UNDEF" in line:
                continue
            if "SECT1" not in line:
                continue
            if "()" not in line:
                continue
            symbol_name = line.split("|")[-1].strip()
            if symbol_name in written_symbols:
                continue
            written_symbols.add(symbol_name)
            yield symbol_name


def symbol_out_to_def_file(symbol_file: pathlib.Path):
    stem_clean = symbol_file.stem.replace(LIB_RAW_PREFIX, "")
    def_file = symbol_file.with_name(f"{stem_clean}.def")
    with open(def_file, "w") as fdef:
        fdef.write("EXPORTS\n")
        for sym_name in iter_symbol_names(symbol_file):
            fdef.write(f"	{sym_name}\n")


def create_symbol_def(module_name: str, def_option: DEFOption):
    symbol_output_file = get_obj_sym_file(module_name)
    if not symbol_output_file.exists():
        create_symbols_file(module_name, def_option)

    symbol_out_to_def_file(symbol_output_file)


if __name__ == "__main__":
    create_symbol_def("bibc", DEFOption.NO_DEF)
    create_symbol_def("bibcxx", DEFOption.NO_DEF)
    create_symbol_def("bibfor", DEFOption.NO_DEF)
    create_symbol_def("aster", DEFOption.NO_DEF)
