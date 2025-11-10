#!/usr/bin/env python3
"""
Generate asterGC.def file for MSVC from compiled object files.

This script extracts symbols from compiled object files in the libs/gc directory
and generates a module definition (.def) file for use with the MSVC linker.

The asterGC library contains both Fortran and C++ code for garbage collection
and memory management.
"""

import argparse
import subprocess
import sys
from pathlib import Path
from def_gen_cache import DumpbinCache, run_dumpbin_for_file


def find_object_files(build_dir):
    """Find all object files in the libs/gc build directory.

    Supports both .obj and .o file extensions.
    """
    obj_files = []
    gc_dir = build_dir / "libs" / "gc"

    if not gc_dir.exists():
        print(f"Warning: libs/gc build directory not found: {gc_dir}")
        return obj_files

    for pattern in ["*.obj", "*.o"]:
        for obj_file in gc_dir.rglob(pattern):
            obj_files.append(obj_file)

    print(f"Found {len(obj_files)} object files in {gc_dir}")
    return obj_files


def extract_symbols_from_dumpbin_lines(lines):
    """Extract and filter symbols from dumpbin output lines.

    Args:
        lines: List of dumpbin output lines

    Returns:
        List of filtered symbol names
    """
    import re

    symbols = []

    # Exclude common non-exportable prefixes and thunks, plus const pools
    exclude_prefixes = (
        "__imp_", "__Cxx", "__RT", "_TI", "_CT", "_Init_thread_", "$", "._", "__chkstk",
        "__real@", "__xmm@", "__ymm@", "__int@", "__m128@", "__m256@", "pybind11"
    )

    # Specific CRT/utility symbols we should not export from AsterGC
    exclude_symbols = {
        "printf", "fprintf", "sprintf", "sprintf_s", "snprintf", "_snprintf",
        "_vfprintf_l", "_vsprintf_l", "_vsprintf_s_l", "_scprintf", "_vscprintf_l",
        "__local_stdio_printf_options",
    }

    # Accept only C-like api names (letters/underscores and digits)
    api_name_re = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")

    for raw in lines:
        line = raw.strip()
        if "External" not in line or "|" not in line:
            continue
        if "UNDEF" in line:
            continue
        if "SECT" not in line:
            continue

        try:
            right = line.split("|", 1)[1].strip()
        except Exception:
            continue
        if not right:
            continue

        token = right.split()[-1]
        name = token.strip()
        if name.endswith("()"):
            name = name[:-2]
        if not name:
            continue

        # Hard rejections
        if name.startswith(exclude_prefixes):
            continue
        # Exclude MSVC and Itanium C++ mangled names
        if name.startswith('?') or name.startswith('_Z'):
            continue
        # Exclude decorated/section names
        if ("@" in name) or ("." in name):
            continue
        if not api_name_re.match(name):
            continue
        if name in exclude_symbols:
            continue
        # Positive allow-list: keep only clear public API prefixes
        allowed_prefixes = ("PyInit_", "aster_gc_", "gc_")
        #if not name.startswith(allowed_prefixes):
        #    continue

        symbols.append(name)

    return symbols


def extract_symbols(obj_files, cache=None, use_hash=False):
    """Extract intended export symbols using dumpbin with robust filtering.

    Goal: avoid exporting MSVC C++ mangled names and template/CRT internals; keep
    only C/Fortran-like API symbols that are actually defined in these objects.

    Args:
        obj_files: List of object file paths
        cache: Optional DumpbinCache instance for caching
        use_hash: If True, use file hash for cache validation; if False, use mtime
    """
    if cache is not None:
        # Use cache-aware extraction
        def extract_for_file(files):
            """Extract symbols from a single file for caching."""
            assert len(files) == 1
            success, lines = run_dumpbin_for_file(files[0])
            if success:
                return extract_symbols_from_dumpbin_lines(lines)
            return []

        return cache.extract_symbols_with_cache(obj_files, extract_for_file, use_hash)
    else:
        # Original non-cached implementation
        symbols = []
        for obj_file in obj_files:
            success, lines = run_dumpbin_for_file(obj_file)
            if success:
                symbols.extend(extract_symbols_from_dumpbin_lines(lines))
        return symbols


def generate_def_file(symbols, output_file):
    """Generate the .def file for asterGC."""
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("LIBRARY AsterGC\n")
        f.write("EXPORTS\n")

        for symbol in symbols:
            f.write(f"    {symbol}\n")

    print(f"Generated {output_file} with {len(symbols)} symbols")


def main():
    parser = argparse.ArgumentParser(description="Generate asterGC.def from object files")
    parser.add_argument("--build-dir", type=Path, help="Build directory containing object files")
    parser.add_argument("--output", type=Path, help="Output .def file path")
    parser.add_argument("--cache", type=Path, help="Cache file path (default: <output-dir>/.astergc_defgen_cache.json)")
    parser.add_argument("--use-hash", action="store_true", help="Use file hash instead of mtime for cache validation (slower but more reliable)")
    parser.add_argument("--no-cache", action="store_true", help="Disable caching")
    args = parser.parse_args()

    # Determine paths
    script_dir = Path(__file__).parent
    project_root = script_dir.parent

    # Use provided build directory or try to find one
    if args.build_dir:
        build_dir = args.build_dir
        print(f"Using provided build directory: {build_dir}")
    else:
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
        print("Error: No object files found. Please compile asterGC first.")
        sys.exit(1)

    # Initialize cache if enabled
    cache = None
    if not args.no_cache:
        if args.cache:
            cache_file = args.cache
        else:
            # Default cache location
            output_file = args.output if args.output else (script_dir / "asterGC.def")
            cache_file = output_file.parent / ".astergc_defgen_cache.json"

        cache = DumpbinCache(cache_file)
        print(f"Using cache file: {cache_file}")

    # Extract symbols
    symbols = sorted(set(extract_symbols(obj_files, cache=cache, use_hash=args.use_hash)))

    # Save cache if used
    if cache is not None:
        cache.save()

    if not symbols:
        print("Warning: No symbols extracted!")
        sys.exit(1)

    # Generate .def file
    output_file = args.output if args.output else (script_dir / "asterGC.def")
    generate_def_file(symbols, output_file)

    print(f"Successfully generated {output_file}")


if __name__ == "__main__":
    main()

