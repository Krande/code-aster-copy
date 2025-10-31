#!/usr/bin/env python3
"""
Generate bibcxx.def file for MSVC from compiled C++ object files.

This script extracts C++ symbols from compiled object files in the bibcxx directory
and generates a module definition (.def) file for use with the MSVC linker.
"""

import subprocess
import sys
from pathlib import Path


def find_object_files(build_dir):
    """Find all object files in the bibcxx build directory.

    Supports both .obj and .o (e.g., .c.1.o) file extensions.
    """
    obj_files = []
    bibcxx_dir = build_dir / "bibcxx"

    if not bibcxx_dir.exists():
        print(f"Warning: bibcxx build directory not found: {bibcxx_dir}")
        return obj_files

    for pattern in ["*.obj", "*.o"]:
        for obj_file in bibcxx_dir.rglob(pattern):
            obj_files.append(obj_file)

    print(f"Found {len(obj_files)} object files in {bibcxx_dir}")
    return obj_files


def extract_symbols(obj_files):
    """Extract C++ symbols from object files using llvm-nm."""
    symbols = set()

    for obj_file in obj_files:
        try:
            # Run llvm-nm to get symbols
            result = subprocess.run(
                ["llvm-nm", "--extern-only", "--defined-only", str(obj_file)],
                capture_output=True,
                text=True,
                check=True
            )

            for line in result.stdout.splitlines():
                parts = line.split()
                if len(parts) < 3:
                    continue

                symbol_type = parts[1]
                symbol_name = parts[2]

                # Filter for C++ symbols
                # C++ mangled names start with '?'
                # T = Text/Code symbol (functions)

                if symbol_type in ['T', 't']:
                    # Include C++ mangled symbols (starting with ?)
                    if symbol_name.startswith('?'):
                        symbols.add(symbol_name)
                    # Also include std:: symbols even if not mangled
                    elif 'std@@' in symbol_name:
                        symbols.add(symbol_name)

        except subprocess.CalledProcessError as e:
            print(f"Warning: Failed to extract symbols from {obj_file}: {e}")
        except FileNotFoundError:
            print("Error: llvm-nm not found. Make sure LLVM tools are in PATH.")
            sys.exit(1)

    return sorted(symbols)


def generate_def_file(symbols, output_file):
    """Generate the .def file."""
    with open(output_file, 'w') as f:
        f.write("LIBRARY bibcxx\n")
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
        print("Error: No object files found. Please compile bibcxx first.")
        sys.exit(1)

    # Extract symbols
    symbols = extract_symbols(obj_files)

    if not symbols:
        print("Warning: No C++ symbols extracted!")
        sys.exit(1)

    # Generate .def file
    output_file = script_dir / "bibcxx.def"
    generate_def_file(symbols, output_file)

    print(f"Successfully generated {output_file}")


if __name__ == "__main__":
    main()

