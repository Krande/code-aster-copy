#!/usr/bin/env python3
"""
Generate bibcxx.def file for MSVC from compiled C++ object files.

This script extracts C++ symbols from compiled object files in the bibcxx directory
and generates a module definition (.def) file for use with the MSVC linker.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path
from def_gen_cache import DumpbinCache, run_dumpbin_for_file


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


def extract_symbols_from_dumpbin_lines(lines):
    """Extract and filter symbols from dumpbin output lines.

    Args:
        lines: List of dumpbin output lines

    Returns:
        List of filtered symbol names
    """
    symbols = []

    # Exclude common non-exportable prefixes and thunks, plus const pools
    exclude_prefixes = (
        "__imp_", "__Cxx", "__RT", "_TI", "_CT", "_Init_thread_", "._", "__chkstk",
        "__real@", "__xmm@", "__ymm@", "__int@", "__m128@", "__m256@",
        # Project/tooling-specific we never want to export from bibcxx
        "pybind11_",
    )

    # Specific CRT/utility symbols we should not export from bibcxx
    exclude_symbols = {
        "printf", "fprintf", "sprintf", "sprintf_s", "snprintf", "_snprintf",
        "_vfprintf_l", "_vsprintf_l", "_vsprintf_s_l", "_scprintf", "_vscprintf_l",
        "__local_stdio_printf_options", "wmemcpy", "wmemset", "CODEASTER_ARRAY_API",
    }

    api_name_re = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")

    # Positive allow-list: restrict bibcxx to known C-API surfaces
    allowed_prefixes = (
        "mgis_",            # MGIS C API
        "vector_vector_",   # small vector wrapper API
        "PyInit_",          # Python module init (if any)
    )
    allowed_exact = {
        "eval_formula_",
        "uexcep_",
        "fma_double_",
        "_raiseException",
        "getMaterialPropertiesNames",  # single legacy helper exported on Linux builds
        # Intentionally do NOT include getTridimMaterialPropertiesNames unless proven defined
    }

    # Explicitly-allowed MSVC C++ mangled symbols that must be exported from bibcxx
    # Keep this list minimal to preserve a stable C-like export surface.
    allowed_msvc_mangled = {
        # std::string toLower(std::string const&)
        "?toLower@@YA?AV?$basic_string@DU?$char_traits@D@std@@V?$allocator@D@2@@std@@AEBV12@@Z",
    }

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

        # Take the FIRST token after '|': this is the real COFF symbol name.
        token = right.split()[0]
        name = token.strip()
        # Drop optional trailing parentheses dumpbin may append for functions
        if name.endswith("()"):
            name = name[:-2]
        if not name:
            continue

        # Check if this is an explicitly allowed mangled symbol first
        if name in allowed_msvc_mangled:
            symbols.append(name)
            continue

        # Hard rejections
        if name in exclude_symbols:
            continue
        if name.startswith(exclude_prefixes):
            continue
        if name.startswith('?') or name.startswith('_Z'):
            # Exclude MSVC and Itanium C++ mangled names (unless in allowed_msvc_mangled)
            continue
        if ("@" in name) or ("." in name):
            continue
        if not api_name_re.match(name):
            continue

        # Positive allow-list enforcement
        if not (name in allowed_exact or any(name.startswith(p) for p in allowed_prefixes)):
            continue

        symbols.append(name)

    return symbols


def extract_symbols(obj_files, cache=None, use_hash=False):
    """Extract intended export symbols using dumpbin with robust filtering.

    Goal: avoid exporting MSVC C++ mangled names and template internals; keep
    only C-like API symbols that are actually defined in these objects.

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


def generate_def_file(symbols, output_file):
    """Generate the .def file."""
    with open(output_file, 'w') as f:
        f.write("LIBRARY bibcxx\n")
        f.write("EXPORTS\n")

        for symbol in symbols:
            f.write(f"    {symbol}\n")

    print(f"Generated {output_file} with {len(symbols)} symbols")


def main():
    parser = argparse.ArgumentParser(description="Generate bibcxx.def from C++ object files")
    parser.add_argument("--build-dir", type=Path, help="Build directory containing object files")
    parser.add_argument("--output", type=Path, help="Output .def file path")
    parser.add_argument("--cache", type=Path, help="Cache file path (default: <output-dir>/.bibcxx_defgen_cache.json)")
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
        print("Error: No object files found. Please compile bibcxx first.")
        sys.exit(1)

    # Initialize cache if enabled
    cache = None
    if not args.no_cache:
        if args.cache:
            cache_file = args.cache
        else:
            # Default cache location
            output_file = args.output if args.output else (script_dir / "bibcxx.def")
            cache_file = output_file.parent / ".bibcxx_defgen_cache.json"

        cache = DumpbinCache(cache_file)
        print(f"Using cache file: {cache_file}")

    # Extract symbols
    symbols = set(extract_symbols(obj_files, cache=cache, use_hash=args.use_hash))

    # Save cache if used
    if cache is not None:
        cache.save()

    # Exclude symbols that belong to separate libs (asterGC, bibfor_ext) if their DEFs exist
    excluded = set()
    excluded |= _read_def_exports(script_dir / "asterGC.def")
    excluded |= _read_def_exports(script_dir / "bibfor_ext.def")
    if excluded:
        before = len(symbols)
        symbols = [s for s in symbols if s not in excluded]
        after = len(symbols)
        print(f"Excluded {before - after} symbols present in separate libs (asterGC/bibfor_ext)")

    if not symbols:
        print("Warning: No symbols left to export after filtering!")
        sys.exit(1)

    symbols = sorted(symbols)

    # Generate .def file
    output_file = args.output if args.output else (script_dir / "bibcxx.def")
    generate_def_file(symbols, output_file)

    print(f"Successfully generated {output_file}")


if __name__ == "__main__":
    main()

