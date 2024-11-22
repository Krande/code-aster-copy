import pathlib

from config import BUILD_DIR, ROOT_DIR, get_bibc_compile_files, TMP_DIR
from msvc_utils import call_using_env
from def_utils import iter_symbol_names


def filter_symbols(def_file: str | pathlib.Path):
    with open(def_file, "r") as fdef:
        export = False
        for line in fdef:
            if line.startswith("EXPORTS"):
                export = True
                continue
            if not export:
                continue
            yield line.strip()


def main():
    TMP_DIR.mkdir(exist_ok=True, parents=True)

    split_files = [
        "supervis/aster_module.c",
        "supervis/aster_core_module.c",
        "supervis/aster_fonctions_module.c",
        "supervis/med_aster_module.c"
    ]
    bibc_o_files = {f.stem.split('.')[0]: f for f in get_bibc_compile_files()}
    original_bibc_def_file = ROOT_DIR / "bibc/bibc.def"
    symbols = list(filter_symbols(original_bibc_def_file))

    for f in split_files:
        fname = f.split("/")[-1].replace(".c", "")
        o_file = bibc_o_files.get(fname)
        symbol_output_file = TMP_DIR / f"{fname}_sym.txt"
        def_name = fname.replace('_module', '') if fname != "aster_module" else fname
        def_file = TMP_DIR / f"{def_name}.def"
        if not symbol_output_file.exists():
            result = call_using_env(["dumpbin", "/symbols", o_file.as_posix(), ">", symbol_output_file], ROOT_DIR)
            if result.returncode != 0:
                raise ValueError(f"Error dumping symbols for {fname} due to {result.stderr=}, {result.stdout=}")
        with open(def_file, "w") as fdef:
            fdef.write(f"LIBRARY {def_name}\n")
            fdef.write("EXPORTS\n")
            for sym_name in iter_symbol_names(symbol_output_file):
                fdef.write(f"	{sym_name}\n")
                if sym_name in symbols:
                    symbols.remove(sym_name)

        new_bibc_def_file = TMP_DIR / "bibc.def"
        with open(new_bibc_def_file, 'w') as fdef:
            fdef.write(f"LIBRARY bibc\n")
            fdef.write("EXPORTS\n")
            for sym_name in symbols:
                if sym_name == '':
                    continue
                fdef.write(f"	{sym_name}\n")


if __name__ == '__main__':
    main()
