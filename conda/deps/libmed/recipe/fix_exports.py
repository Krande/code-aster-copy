#!/usr/bin/env python3
"""Regenerate import libraries with lowercase+underscore aliases.

When libmed is compiled with ifx defaults (UPPERCASE Fortran symbols), the
import libraries (.lib) export e.g. 'MFACRE'. But code_aster, compiled with
/names:lowercase /assume:underscore, expects 'mfacre_'.

This script:
1. Reads exported symbols from installed DLLs using dumpbin /exports
2. Creates a .def file with both original exports and lowercase_ aliases
3. Regenerates the .lib import library using lib.exe /DEF
"""

import os
import re
import subprocess
import sys


def get_dll_exports(dll_path):
    """Get exported symbol names from a DLL using dumpbin."""
    result = subprocess.run(
        ["dumpbin", "/exports", dll_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"  Warning: dumpbin failed for {dll_path}: {result.stderr}")
        return []

    symbols = []
    # dumpbin /exports output has a table with columns:
    # ordinal  hint  RVA  name
    # We look for lines with 4+ columns where the 4th is the symbol name
    in_table = False
    for line in result.stdout.splitlines():
        stripped = line.strip()
        if stripped.startswith("ordinal"):
            in_table = True
            continue
        if in_table and stripped == "":
            # End of table (blank line after entries)
            if symbols:
                break
            continue
        if in_table:
            parts = stripped.split()
            if len(parts) >= 4:
                # ordinal hint RVA name
                sym = parts[3]
                symbols.append(sym)
            elif len(parts) == 2:
                # Some entries have: ordinal name (forwarded)
                sym = parts[1]
                if not sym.startswith("("):
                    symbols.append(sym)
    return symbols


def create_def_with_aliases(dll_name, symbols):
    """Create a .def file content with original exports and lowercase_ aliases."""
    lines = [f"LIBRARY {dll_name}", "EXPORTS"]
    existing = set(symbols)

    for sym in symbols:
        lines.append(f"    {sym}")
        alias = sym.lower() + "_"
        if alias != sym and alias not in existing:
            lines.append(f"    {alias} = {sym}")
            existing.add(alias)

    return "\n".join(lines) + "\n"


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <bin_dir> <lib_dir>")
        print("  bin_dir: directory containing installed DLLs")
        print("  lib_dir: directory containing installed .lib files")
        sys.exit(1)

    bin_dir = sys.argv[1]
    lib_dir = sys.argv[2]

    # Target libraries that need Fortran aliases
    targets = ["medfwrap", "med"]

    total_aliases = 0
    for target in targets:
        dll_path = os.path.join(bin_dir, f"{target}.dll")
        lib_path = os.path.join(lib_dir, f"{target}.lib")

        if not os.path.exists(dll_path):
            # Try with 'lib' prefix
            dll_path = os.path.join(bin_dir, f"lib{target}.dll")
            lib_path = os.path.join(lib_dir, f"lib{target}.lib")

        if not os.path.exists(dll_path):
            print(f"  Skipping {target}: DLL not found")
            continue

        print(f"  Processing {dll_path}")
        symbols = get_dll_exports(dll_path)
        if not symbols:
            print(f"  Warning: no exports found in {dll_path}")
            continue

        dll_name = os.path.basename(dll_path)
        def_content = create_def_with_aliases(dll_name, symbols)

        # Write temporary .def file
        def_path = os.path.join(lib_dir, f"{target}_aliases.def")
        with open(def_path, "w") as f:
            f.write(def_content)

        n_aliases = sum(1 for s in symbols if s.lower() + "_" != s)
        total_aliases += n_aliases
        print(f"  {target}: {len(symbols)} exports, {n_aliases} aliases added")

        # Regenerate .lib using lib.exe /DEF
        result = subprocess.run(
            ["lib.exe", f"/DEF:{def_path}", f"/OUT:{lib_path}", "/MACHINE:X64"],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            print(f"  Error: lib.exe failed for {target}: {result.stderr}")
            sys.exit(1)

        print(f"  Regenerated {lib_path}")

    print(f"Total: {total_aliases} aliases across {len(targets)} libraries")


if __name__ == "__main__":
    main()
