#!/usr/bin/env python3
"""
Generate bibfor.def file for MSVC from compiled Fortran object files.

This script extracts Fortran symbols from compiled object files in the bibfor directory
and generates a module definition (.def) file for use with the MSVC linker.
"""

import subprocess
import sys
from pathlib import Path


def find_object_files(build_dir):
    """Find all object files in the bibfor build directory.

    Supports both .obj and .o (e.g., .f.1.o) file extensions.
    """
    obj_files = []
    bibfor_dir = build_dir / "bibfor"

    if not bibfor_dir.exists():
        print(f"Warning: bibfor build directory not found: {bibfor_dir}")
        return obj_files

    for pattern in ["*.obj", "*.o"]:
        for obj_file in bibfor_dir.rglob(pattern):
            obj_files.append(obj_file)

    print(f"Found {len(obj_files)} object files in {bibfor_dir}")
    return obj_files


def extract_symbols(obj_files):
    """Extract defined external symbols from COFF objs via dumpbin.

    We no longer require trailing underscores nor the "()" marker, because
    MSVC/Intel Fortran often emit symbols without underscores and may not
    annotate with parentheses in COFF symbol listings.

    Additionally, filter out MSVC/Intel Fortran constant pool symbols like
    __real@..., __xmm@..., __ymm@..., and similar, which must not be exported.
    """
    import re

    defined = set()

    # Exclude common non-exportable prefixes and thunks
    exclude_prefixes = (
        "__imp_", "__Cxx", "__RT", "_TI", "_CT", "_Init_thread_", "$", "._", "__chkstk",
        # Constant pool / vectorized literals and others we must not export
        "__real@", "__xmm@", "__ymm@", "__int@", "__m128@", "__m256@",
    )

    # Accept only reasonable API-like symbols: start with a letter, then [A-Za-z0-9_]*
    # This will drop names containing '@' (stdcall decorations), dots, etc.
    api_name_re = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")

    for obj_file in obj_files:
        try:
            result = subprocess.run(
                ["dumpbin", "/SYMBOLS", str(obj_file)],
                capture_output=True,
                text=True,
                check=True,
            )
        except FileNotFoundError:
            print("Error: dumpbin.exe not found. Make sure MSVC tools are in PATH.")
            sys.exit(1)
        except subprocess.CalledProcessError:
            continue

        for raw in result.stdout.splitlines():
            line = raw.strip()
            # We want: External, defined (not UNDEF), associated with a section
            if "External" not in line or "|" not in line:
                continue
            if "UNDEF" in line:
                continue
            if "SECT" not in line:
                continue

            # right side after the '|'
            try:
                right = line.split("|", 1)[1].strip()
            except Exception:
                continue
            if not right:
                continue

            # The symbol token is typically the last field after the '|'
            token = right.split()[-1]
            name = token.strip()
            # Some dumpbin variants append '()' to function names; strip it if present
            if name.endswith("()"):
                name = name[:-2]
            if not name:
                continue

            # Exclusions for obvious non-APIs and internals
            if name.startswith(exclude_prefixes):
                continue
            if name.startswith("?"):  # C++ mangled
                continue
            if "." in name:          # module metadata
                continue
            if name.endswith("._") or "_module" in name or "_mp_" in name:
                continue
            if "@" in name:          # stdcall-size decorated or const pools
                continue
            if not api_name_re.match(name):
                continue

            defined.add(name)

    return sorted(defined)


def generate_def_file(symbols, output_file):
    """Generate a .def exporting both MSVC/ifx and GNU-style names via aliases.

    Rules per discovered internal symbol `s`:
    - Always export the internal symbol `s` as-is.
    - If `s` ends with `_` (GNU style), also export the alias without underscore: `base = s[:-1]`; emit `base=s`.
    - If `s` does NOT end with `_` (MSVC/ifx style), also export the alias with underscore: emit `s_=s`.

    This ensures both spellings are available regardless of which one is actually defined internally.
    """
    exports = []
    for s in symbols:
        exports.append(s)  # always export the actual internal symbol
        if s.endswith("_"):
            base = s[:-1]
            if base:
                exports.append(f"{base}={s}")
        else:
            exports.append(f"{s}_={s}")

    # de-duplicate while preserving order
    seen = set()
    ordered = []
    for e in exports:
        if e not in seen:
            seen.add(e)
            ordered.append(e)

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("LIBRARY bibfor\n")
        f.write("EXPORTS\n")
        for e in ordered:
            f.write(f"    {e}\n")

    print(f"Generated {output_file} with {len(ordered)} exports (from {len(symbols)} symbols)")


def main():
    # Determine paths
    script_dir = Path(__file__).parent
    project_root = script_dir.parent

    # Try to find build directory
    build_dirs = []
    for pattern in ["build/int64/debug", "build/int64/release", "build/int32/debug", "build/int32/release"]:
        build_path = project_root / pattern
        if build_path.exists():
            build_dirs.append(build_path)

    if not build_dirs:
        print("Error: No build directory found. Please run 'waf build' first.")
        sys.exit(1)

    # Use the first found build directory (most recent)
    build_dir = build_dirs[0]
    print(f"Using build directory: {build_dir}")

    # Find object files
    obj_files = find_object_files(build_dir)

    if not obj_files:
        print("Error: No object files found. Please compile bibfor first.")
        sys.exit(1)

    # Extract symbols
    symbols = extract_symbols(obj_files)

    if not symbols:
        print("Warning: No Fortran symbols extracted!")
        sys.exit(1)

    # Generate .def file
    output_file = script_dir / "bibfor.def"
    generate_def_file(symbols, output_file)

    print(f"Successfully generated {output_file}")


if __name__ == "__main__":
    main()

