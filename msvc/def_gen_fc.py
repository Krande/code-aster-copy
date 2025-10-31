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
    """Extract Fortran symbols from object files using dumpbin (COFF-aware).

    Rules:
    - Only External symbols that are defined in this object (exclude UNDEF)
    - Exclude C++/CRT/RTTI/import thunks and similar internals
    - Keep only names ending with '_' (Fortran convention)
    - Exclude any names containing '.' (module metadata like module._)
    - Exclude obvious module-related exports ending with '._' or containing '_module'
    """
    symbols = set()

    exclude_prefixes = (
        "__imp_", "__Cxx", "__RT", "_TI", "_CT", "_Init_thread_", "$", "._", "__chkstk"
    )

    for obj_file in obj_files:
        try:
            result = subprocess.run(
                ["dumpbin", "/SYMBOLS", str(obj_file)],
                capture_output=True,
                text=True,
                check=True
            )
            for line in result.stdout.splitlines():
                if "External" not in line or "|" not in line:
                    continue
                if "UNDEF" in line:
                    continue
                # Only include entries that are associated with a real section and look like functions
                if "SECT" not in line:
                    continue
                if "()" not in line:
                    continue
                parts = line.split('|', 1)
                if len(parts) < 2:
                    continue
                right = parts[1].strip()
                if not right:
                    continue
                token = right.split()[0]
                name = token.strip().strip('()')

                if not name:
                    continue
                if name.startswith(exclude_prefixes):
                    continue
                if name.startswith('?'):
                    continue
                if not name.endswith('_'):
                    continue
                if '.' in name:
                    continue
                if name.endswith('._') or '_module' in name or "_mp_" in name:
                    continue

                symbols.add(name)
        except subprocess.CalledProcessError:
            continue
        except FileNotFoundError:
            print("Error: dumpbin.exe not found. Make sure MSVC tools are in PATH.")
            sys.exit(1)

    return sorted(symbols)


def generate_def_file(symbols, output_file):
    """Generate the .def file."""
    with open(output_file, 'w') as f:
        f.write("LIBRARY bibfor\n")
        f.write("EXPORTS\n")

        for symbol in symbols:
            f.write(f"    {symbol}\n")

    print(f"Generated {output_file} with {len(symbols)} symbols")


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

