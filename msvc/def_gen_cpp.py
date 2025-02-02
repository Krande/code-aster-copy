#!/usr/bin/env python3

import argparse
import os
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
    def __init__(self):
        prefix = os.getenv("CONDA_PREFIX")

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

        self.clang_cl = clang_cl_exe
        self.llvm_nm = llvm_nm_exe
        self.prefix = prefix
        bibc_include_paths = ROOT_DIR / "bibc"
        bibcxx_include_paths  = ROOT_DIR / "bibcxx"
        py_include = prefix / "include"
        self.bibcxx_include_paths = [f"/I{bibcxx_include_paths.as_posix()}", f"/I{bibcxx_include_paths.as_posix()}/include", f"/I{bibc_include_paths.as_posix()}/include", f"/I{py_include.as_posix()}"]

        self.symbols_filter = [
            "Ref_count_obj",
            "shared_ptr@",
            "Incref@",
            "Destroy@",
            "Delete_this@",
            "Decref@"
        ]
        # make symbols lower
        self.symbols_filter = [f.lower() for f in self.symbols_filter]

    def should_skip_symbol(self, symbol: str) -> bool:
        result = any(f in symbol.lower() for f in self.symbols_filter)
        return result


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
            if scanner.should_skip_symbol(name):
                continue
            if sym_type.isupper() and sym_type != "U":
                symbols.add(name)
    return symbols

def build_def(sources,obj_dir,lib,output, use_as_cli=False, compile_objs=False) -> None:
    if use_as_cli:
        args = parse_args()
        sources = args.sources
        obj_dir = args.obj_dir
        lib = args.lib
        output = args.output

    scanner = Scanner()
    all_symbols: set[str] = set()
    if not compile_objs:
        obj_files = Path("../build/std/debug/bibcxx").resolve().absolute().rglob("*.o")
        for obj_file in obj_files:
            symbols = extract_symbols(obj_file, scanner)
            all_symbols |= symbols
    else:
        for src in sources:
            obj_file = compile_source(src, obj_dir, scanner)
            symbols = extract_symbols(obj_file, scanner)
            all_symbols |= symbols

    print(f"Found {len(all_symbols)} symbols.")

    # Write the .def file with the standard format.
    with output.open("w") as f:
        f.write(f"LIBRARY {lib}\n")
        f.write("EXPORTS\n")
        for sym in sorted(all_symbols):
            f.write(f"    {sym}\n")

    print(f".def file written to {output}")


def main():
    sources = Path("../bibcxx").resolve().absolute().rglob("*.cxx")
    obj_dir = Path("temp")
    lib = "bibcxx"
    output = Path("bibcxx.def")
    obj_dir.mkdir(parents=True, exist_ok=True)

    build_def(sources, obj_dir,lib,output,False)


if __name__ == "__main__":
    main()
