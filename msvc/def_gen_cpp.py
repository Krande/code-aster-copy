#!/usr/bin/env python3

import argparse
import json
import os
import pathlib
import shutil
import subprocess
from pathlib import Path

ROOT_DIR = Path(__file__).parent.parent

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a .def file with mangled symbols from C++ sources using clangdev"
    )
    parser.add_argument(
        "sources",
        nargs="+",
        type=Path,
        help="C++ source files to process"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("exports.def"),
        help="Output .def file (default: exports.def)"
    )
    parser.add_argument(
        "--lib",
        type=str,
        default="generated",
        help="Library name to use in the .def file (default: generated)"
    )
    parser.add_argument(
        "--obj-dir",
        type=Path,
        default=Path("obj"),
        help="Directory to store intermediate object files (default: obj)"
    )
    args = parser.parse_args()
    args.obj_dir.mkdir(parents=True, exist_ok=True)
    return args


class Scanner:
    def __init__(self, bib_name, bib_dir):
        prefix = os.getenv("CONDA_PREFIX")
        self.bib_name = bib_name
        if not prefix:  # Probably running from IDE
            prefix = (ROOT_DIR / ".pixi/envs/dev").resolve().absolute()
            llvm_nm_exe = shutil.which("llvm-nm", path=prefix / "Library/bin")
            clang_cl_exe = shutil.which("clang-cl", path=prefix / "Library/bin")
        else:  # all envs are set as they should and llvm-nm should be in the path
            llvm_nm_exe = shutil.which("llvm-nm")
            clang_cl_exe = shutil.which("clang-cl")

        if not llvm_nm_exe:
            raise RuntimeError("llvm-nm not found in PATH")

        if not clang_cl_exe:
            raise RuntimeError("clang-cl not found in PATH")

        self.bib_dir = Path(bib_dir)
        if not self.bib_dir.exists():
            raise RuntimeError(f"{bib_name} not found in build directory {bib_dir}")

        self.clang_cl = clang_cl_exe
        self.llvm_nm = llvm_nm_exe
        self.prefix = prefix

        bibc_include_paths = ROOT_DIR / "bibc"
        bibcxx_include_paths  = ROOT_DIR / "bibcxx"
        py_include = prefix / "include"
        self.bibcxx_include_paths = [f"/I{bibcxx_include_paths.as_posix()}", f"/I{bibcxx_include_paths.as_posix()}/include", f"/I{bibc_include_paths.as_posix()}/include", f"/I{py_include.as_posix()}"]

        self.symbols_containing_string_filter = [
            "Ref_count",
            "shared_ptr@",
            "Incref@",
            "Destroy@",
            "Delete_this@",
            "Decref@",
            "_vfprintf_l",
            "_vscprintf_l",
            "wmemcpy",
            "wmemset",
            "bad_alloc",
            "bad_cast",
            "?use_count@?$_Ptr_base@",  # Internal use_count function.
            "ErrorCpp",
            "pybind11@@" # Internal pybind11 symbols.

        ]
        # make symbols lower
        self.symbols_containing_string_filter = [f.lower() for f in self.symbols_containing_string_filter]
        self.symbols_starting_with_filter = [
            "??_C@", # Internal C++ runtime symbols.
            "??_G",  # Internal C++ runtime symbols.
            "_CT",  # Internal C++ runtime symbols.
            "_TI",  # Internal C++ runtime symbols.
            "__",
            "??_8"
        ]
        self.symbols_starting_with_filter = [f.lower() for f in self.symbols_starting_with_filter]

    def should_skip_symbol_based_on_contained_string(self, symbol: str) -> bool:
        result = any(f in symbol.lower() for f in self.symbols_containing_string_filter)
        return result

    def should_skip_symbol_based_on_starting_string(self, symbol: str) -> bool:
        result = any(symbol.lower().startswith(f) for f in self.symbols_starting_with_filter)
        return result

    def compile_sources(self, sources, obj_dir):
        for src in sources:
            obj_file = compile_source(src, obj_dir, self)

def compile_source(source: Path, obj_dir: Path, scanner: Scanner) -> Path:
    """Compile a C++ source file to an object file using clangdev."""
    obj_path = obj_dir / (source.stem + ".obj")
    print(f"Compiling {source} -> {obj_path}")
    subprocess.run(
        [scanner.clang_cl,*scanner.bibcxx_include_paths, "-c", str(source), "-o", str(obj_path)],
        check=True
    )
    return obj_path

def extract_symbols(obj_file: Path, scanner: Scanner) -> set[str]:
    """
    Run llvm-nm on the object file and return a set of defined mangled symbols.

    This function assumes llvm-nm prints lines in the form:
        <address> <symbol_type> <symbol_name>
    and that symbols with a uppercase symbol type (other than 'U' for undefined)
    are the ones to be exported.
    """
    symbols: set[str] = set()
    print(f"Extracting symbols from {obj_file}")

    result = subprocess.run(
        [scanner.llvm_nm, str(obj_file)],
        capture_output=True,
        text=True,
        shell=True
    )
    if result.returncode != 0:
        print(result.stderr.decode())
        raise RuntimeError("llvm-nm failed")

    for line in result.stdout.splitlines():
        parts = line.strip().split()
        if len(parts) == 3:
            _address, sym_type, name = parts
            if sym_type.isupper() and sym_type != "U":
                symbols.add(name)
    return symbols

def build_def(scanner: Scanner, output, clear_cache=False) -> None:
    all_symbols: set[str] = set()
    cache_file = pathlib.Path("temp/symbols_cache.json")
    if cache_file.exists() and clear_cache is False:
        using_cache = True
        with cache_file.open("r") as f:
            cache_json = json.load(f)
    else:
        cache_json = {}
        using_cache = False

    obj_files = scanner.bib_dir.rglob("*.o")

    for obj_file in obj_files:
        if scanner.bib_name == "bibcxx":
            skip_subdir = scanner.bib_dir / "Mfront"
            if obj_file.parent == skip_subdir:
                continue
        rel_obj_fp = obj_file.relative_to(ROOT_DIR).as_posix()
        if cache_file.exists() and rel_obj_fp in cache_json:
            symbols = cache_json[rel_obj_fp]
        else:
            symbols = extract_symbols(obj_file, scanner)
        cache_json[rel_obj_fp] = []
        filtered_symbols = set()
        for name in symbols:
            if using_cache is False:
                cache_json[rel_obj_fp].append(name)
            if scanner.should_skip_symbol_based_on_contained_string(name):
                continue
            if scanner.should_skip_symbol_based_on_starting_string(name):
                continue
            filtered_symbols.add(name)
        all_symbols |= filtered_symbols


    print(f"Found {len(all_symbols)} symbols.")
    if len(list(all_symbols)) == 0:
        raise RuntimeError("No symbols files found. Please compile the sources first.")

    # Write the .def file with the standard format.
    with output.open("w") as f:
        f.write(f"LIBRARY {scanner.bib_name}\n")
        f.write("EXPORTS\n")
        for sym in sorted(all_symbols):
            f.write(f"    {sym}\n")

    print(f".def file written to {output}")

    if using_cache is False:
        with cache_file.open("w") as f:
            json.dump(cache_json, f)

def main():
    bib_dir = (ROOT_DIR / "build/int64/debug/bibcxx").resolve().absolute()
    scanner = Scanner("bibcxx", bib_dir)

    build_def(scanner, Path("bibcxx.def"), clear_cache=True)


if __name__ == "__main__":
    main()
