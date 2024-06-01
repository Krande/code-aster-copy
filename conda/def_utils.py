import pathlib


def iter_symbol_names(symbol_file: pathlib.Path | str):
    # Based on postgres https://github.com/postgres/postgres/blob/master/src/tools/msvc_gendef.pl
    written_symbols = set()
    header_passed = False
    with open(symbol_file) as f:
        for line in f:
            if "UNDEF" in line:
                continue
            if "COFF SYMBOL TABLE" in line:
                header_passed = True
                continue
            if not header_passed:
                continue
            # if "SECT1" not in line:
            #     continue
            # if "()" not in line:
            #     continue
            if 'Static' in line:
                continue
            if 'Filename' in line:
                continue
            if '|' not in line:
                continue

            symbol_name = line.split("|")[-1].strip()
            if "(`string')" in symbol_name:
                continue
            if symbol_name in written_symbols:
                continue
            if symbol_name.startswith("@"):
                continue
            if symbol_name.startswith(r"\("):
                continue
            if symbol_name.startswith("__real"):
                continue
            if symbol_name.startswith("__xmm"):
                continue
            if symbol_name.startswith("__imp"):
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
