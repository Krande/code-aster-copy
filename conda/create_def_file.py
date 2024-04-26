import pathlib
import subprocess

THIS_DIR = pathlib.Path(__file__).resolve().parent
ROOT_DIR = THIS_DIR.parent
BUILD_DEBUG_DIR = ROOT_DIR / "build/std/debug"
BAT_CALL_FILE = ROOT_DIR / "call_compile.bat"
PREFIX = "tmp_"


def generate_lib_file(lib_name: str) -> pathlib.Path:
    bib_debug_dir = BUILD_DEBUG_DIR / lib_name
    bib_lib_in_file = bib_debug_dir / f"{PREFIX}{lib_name}_lib.in"
    bib_symbols_file = THIS_DIR / f"{PREFIX}{lib_name}.lib"

    with open(bib_lib_in_file, "w") as f:
        f.write(f"/OUT:{bib_symbols_file}\n")
        for fp in bib_debug_dir.rglob("*.o"):
            f.write(f"{fp}\n")

    command = f"call {BAT_CALL_FILE} lib.exe @{bib_lib_in_file}"

    # Execute the command
    result = subprocess.run(command, cwd=ROOT_DIR, shell=True, capture_output=True, text=True)

    # Print output and error for debugging
    print("stdout:", result.stdout)
    print("stderr:", result.stderr)
    return bib_symbols_file


def create_symbols_file(module_name: str) -> None:
    symbol_output_file = (THIS_DIR / f"{PREFIX}{module_name}").with_suffix(".txt")
    bib_lib_file = symbol_output_file.with_suffix(".lib")
    if not bib_lib_file.exists():
        bib_lib_file = generate_lib_file(module_name)
    command = f"call {BAT_CALL_FILE} dumpbin /symbols {bib_lib_file} > {symbol_output_file}"
    result = subprocess.run(command, cwd=ROOT_DIR, shell=True, capture_output=True, text=True)
    print("stdout:", result.stdout)
    print("stderr:", result.stderr)


def symbol_out_to_def_file(symbol_file: pathlib.Path):
    stem_clean = symbol_file.stem.replace(PREFIX, '')
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
                if '00000000' not in line:
                    continue
                symbol_name = line.split("|")[-1].strip()
                if symbol_name in written_symbols:
                    continue
                written_symbols.add(symbol_name)
                fdef.write(f"{symbol_name}\n")


def create_symbol_def(module_name: str):
    symbol_output_file = (THIS_DIR / f"{PREFIX}{module_name}").with_suffix(".txt")
    if not symbol_output_file.exists():
        create_symbols_file(module_name)
    symbol_out_to_def_file(symbol_output_file)


if __name__ == "__main__":
    for mod in ["bibc", "bibfor", "bibcxx"]:
        create_symbol_def(mod)
