#!/usr/bin/env python3
"""
Generate DEF file for bibc (C library) by extracting symbols from compiled object files.

This script uses dumpbin.exe to extract symbols from C object files and creates
a .def file for exporting them from the DLL.
"""

import argparse
import subprocess
import sys
import re
from pathlib import Path


def find_object_files(build_dir, library_name="bibc"):
    """Find all object files for the given library.

    Restrict search to the specific library subdirectory to avoid cross-picking
    from similarly-named directories (e.g., bibcxx). Supports both .obj and .o.
    """
    obj_files = []
    build_path = Path(build_dir)
    lib_dir = build_path / library_name

    if not lib_dir.exists():
        # Fall back to whole tree scan as a last resort (older layouts), but be strict on path
        search_roots = [lib_dir, build_path]
    else:
        search_roots = [lib_dir]

    patterns = ["*.obj", "*.o"]
    for root in search_roots:
        if not root.exists():
            continue
        for pattern in patterns:
            for obj_file in root.rglob(pattern):
                # Ensure path actually contains the library dir name as a segment
                if library_name in obj_file.parts:
                    obj_files.append(obj_file)

    # Deduplicate while preserving order
    seen = set()
    unique_files = []
    for p in obj_files:
        if p not in seen:
            seen.add(p)
            unique_files.append(p)

    return unique_files


def extract_c_symbols(obj_file):
    """Extract C symbols from an object file using dumpbin.

    Notes:
    - Filters out import thunks (__imp_), C++/RTTI/CRT internals, and stray annotations
      appended by dumpbin like "(__declspec(dllimport) ... )".
    """
    symbols = []
    data_symbols = []

    # Prefixes to exclude from DEF (CRT/RTTI/imports/etc.)
    exclude_prefixes = (
        "__imp_", "__Cxx", "__RT", "_TI", "_CT", "_Init_thread_", "$", "._", "__chkstk",
        # Constant pool / vectorized literals and others we must not export
        "__real@", "__xmm@", "__ymm@", "__int@", "__m128@", "__m256@",
    )
    # Specific symbols to exclude (not to be exported from bibc)
    exclude_symbols = {
        "_unused_main_", "_main", "main", "WinMain", "_WinMain@16",
        "DllMain", "_DllMain@12", "DllMainCRTStartup", "_DllMainCRTStartup@12",
        # CRT internal functions that should not be exported
        "_vfprintf_l", "_vsnprintf", "_vsnprintf_l", "_vsprintf_l",
        "_fprintf_l", "_sprintf_l", "_printf_l", "_snprintf_l",
        "_scanf_l", "_fscanf_l", "_sscanf_l",
    }

    try:
        # Run dumpbin to get symbols
        result = subprocess.run(
            ["dumpbin", "/SYMBOLS", str(obj_file)],
            capture_output=True,
            text=True,
            check=True
        )

        for line in result.stdout.split('\n'):
            # Look for External symbols that are actually defined in this object
            if "External" in line and "|" in line:
                if "UNDEF" in line:
                    # skip undefined externals; they are imported from other libs
                    continue
                parts = line.split('|', 1)
                if len(parts) >= 2:
                    # Right side typically starts with the symbol then whitespace and annotations
                    right = parts[1].strip()
                    if not right:
                        continue
                    token = right.split()[0]  # strip any trailing annotations

                    # strip any surrounding parentheses remnants
                    symbol = token.strip() \
                        .removesuffix(')') \
                        .removeprefix('(')

                    # Check if it's a DATA symbol (global variable)
                    is_data = "SECT" in line and not ("()" in line or "notype" in line.lower())

                    # Filtering rules
                    if not symbol:
                        continue
                    if symbol in exclude_symbols:
                        continue
                    # Exclude any symbol with obvious decoration/metadata
                    if symbol.startswith('.') or ('@' in symbol) or ('.' in symbol):
                        continue
                    if symbol.startswith(exclude_prefixes):
                        continue

                    # Exclude C++ mangled names (MSVC and Itanium)
                    if symbol.startswith('?') or symbol.startswith('_Z'):
                        continue

                    # Accept standard C/Fortran-identifiers (letters/underscore, digits/underscore)
                    # This allows lowercase exports like r8prem_ used by Fortran-callable C shims.
                    if not re.match(r'^_?[A-Za-z][A-Za-z0-9_]*$', symbol):
                        continue

                    if is_data:
                        data_symbols.append(symbol)
                    else:
                        symbols.append(symbol)

    except subprocess.CalledProcessError as e:
        print(f"Warning: Failed to process {obj_file}: {e}", file=sys.stderr)
    except FileNotFoundError:
        print("Error: dumpbin.exe not found. Ensure MSVC tools are in PATH.", file=sys.stderr)
        sys.exit(1)

    return symbols, data_symbols


def _read_def_exports(def_path: Path):
    """Read a .def file and return a set of exported symbol names (first token per line after EXPORTS)."""
    exports = set()
    p = Path(def_path)
    if not p.exists():
        return exports
    try:
        lines = p.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return exports
    in_exports = False
    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        up = line.upper()
        if up.startswith("EXPORTS"):
            in_exports = True
            continue
        if not in_exports:
            continue
        token = line.split()[0]
        if token and token.upper() not in ("EXPORTS", "LIBRARY"):
            exports.add(token)
    return exports


def generate_def_file(symbols, data_symbols, output_file, library_name="bibc"):
    """Generate a .def file from the list of symbols."""
    # Remove duplicates and sort
    unique_symbols = sorted(set(symbols))
    unique_data_symbols = sorted(set(data_symbols))

    with open(output_file, 'w') as f:
        f.write(f"LIBRARY {library_name}\n")
        f.write("EXPORTS\n")

        # Write DATA symbols first
        for symbol in unique_data_symbols:
            f.write(f"\t{symbol} \t DATA\n")

        # Write regular symbols
        for symbol in unique_symbols:
            f.write(f"\t{symbol}\n")

    total = len(unique_symbols) + len(unique_data_symbols)
    print(f"Generated {output_file} with {total} symbols ({len(unique_data_symbols)} DATA)")


def main():
    parser = argparse.ArgumentParser(
        description="Generate DEF file for bibc C library"
    )
    parser.add_argument(
        "--build-dir",
        default="build/int64",
        help="Build directory containing object files"
    )
    parser.add_argument(
        "--output",
        default="msvc/bibc.def",
        help="Output DEF file path"
    )
    parser.add_argument(
        "--obj-files",
        nargs='+',
        help="Specific object files to process (optional)"
    )

    args = parser.parse_args()

    # Get object files
    if args.obj_files:
        obj_files = [Path(f) for f in args.obj_files]
    else:
        obj_files = find_object_files(args.build_dir, "bibc")

    if not obj_files:
        print(f"No object files found in {args.build_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(obj_files)} object files...")

    # Extract symbols from all object files
    all_symbols = []
    all_data_symbols = []
    for obj_file in obj_files:
        symbols, data_symbols = extract_c_symbols(obj_file)
        all_symbols.extend(symbols)
        all_data_symbols.extend(data_symbols)

    # Exclude symbols that belong to separate libs (asterGC, bibfor_ext) if their DEFs exist
    script_dir = Path(__file__).parent
    excluded = set()
    excluded |= _read_def_exports(script_dir / "asterGC.def")
    excluded |= _read_def_exports(script_dir / "bibfor_ext.def")
    if excluded:
        before = (len(set(all_symbols)) + len(set(all_data_symbols)))
        all_symbols = [s for s in all_symbols if s not in excluded]
        all_data_symbols = [s for s in all_data_symbols if s not in excluded]
        after = (len(set(all_symbols)) + len(set(all_data_symbols)))
        print(f"Excluded {before - after} symbols present in separate libs (asterGC/bibfor_ext)")

    # Generate DEF file
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    generate_def_file(all_symbols, all_data_symbols, output_path, "bibc")


if __name__ == "__main__":
    main()

