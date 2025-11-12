# DEF File Auto-Generation Integration Plan

## Overview

This document outlines the integration of automatic `.def` file generation into the build process for the MSVC port of code_aster. The `.def` files are required for exporting symbols from shared libraries on Windows using the MSVC toolchain.

## Background

On Windows with MSVC, symbols must be explicitly exported from DLLs using `.def` (module definition) files. These files list all the symbols that should be exported from the shared library.

## Libraries Requiring .def Files

1. **bibfor.def** - Fortran library (contains Fortran subroutines with trailing underscores)
2. **bibcxx.def** - C++ library (contains C++ mangled names)
3. **bibc.def** - C library (contains C functions and Python module init functions)

## Integration Strategy

### Phase 1: Object File Analysis

Before linking each library, the build system will:

1. Collect all compiled object files (.obj) for the library
2. Use `llvm-nm` or `dumpbin` to extract exported symbols from object files
3. Filter symbols based on library type (Fortran vs C vs C++)
4. Generate the `.def` file

### Phase 2: Linking with .def Files

The generated `.def` files are then passed to the linker when building the DLLs.

## Implementation Components

### 1. Symbol Extraction Scripts

#### `msvc/def_gen_fc.py`
- Extracts Fortran symbols (functions ending with `_`)
- Uses `llvm-nm` to list symbols from `.obj` files
- Filters for TEXT (T) symbols that are public
- Generates `bibfor.def`

#### `msvc/def_gen_cpp.py`
- Extracts C++ symbols (mangled names)
- Uses `llvm-nm` to list symbols
- Filters C++ mangled names (starting with `?`)
- Generates `bibcxx.def`

#### `msvc/def_gen_c.py`
- Extracts C symbols
- Filters for non-mangled, non-Fortran symbols
- Includes Python module init functions (`PyInit_*`)
- Generates `bibc.def`

### 2. WAF Integration

The `.def` generation is integrated into the WAF build system:

- Add a pre-link task that runs before `stlib` or `shlib` tasks
- The task scans compiled object files
- Generates the appropriate `.def` file
- The linker task uses the generated `.def` file

### 3. Pixi Tasks

Manual generation tasks are available for debugging:

```cmd
pixi run def-bibfor     # Generate bibfor.def
pixi run def-bibcpp     # Generate bibcxx.def
pixi run def-bibc       # Generate bibc.def
```

## Build Flow

```
Compile Sources (.f, .cpp, .c)
        ↓
Generate Object Files (.obj)
        ↓
[PRE-LINK] Extract Symbols → Generate .def file
        ↓
Link with .def file → Generate DLL
```

## Symbol Filtering Rules

### Fortran (bibfor.def)
- Include: Symbols ending with `_` (Fortran name mangling)
- Exclude: C++ mangled symbols (starting with `?`)
- Exclude: Compiler-generated symbols

### C++ (bibcxx.def)
- Include: C++ mangled symbols (starting with `?`)
- Exclude: Compiler-generated internal symbols
- Exclude: Weak symbols

### C (bibc.def)
- Include: Regular C symbols
- Include: `PyInit_*` functions (Python module initialization)
- Include: DATA symbols (global variables used by Python bindings)
- Exclude: Fortran symbols (ending with `_`)
- Exclude: C++ mangled symbols

## Tools Used

- **llvm-nm**: Symbol extraction from object files (available in clang toolchain)
- **PowerShell**: Script execution on Windows
- **Python**: Symbol filtering and `.def` file generation

## Future Improvements

1. Incremental .def generation (only regenerate when object files change)
2. Symbol visibility annotations in source code
3. Automatic detection of new symbols
4. Cross-check with Linux symbol exports

## Notes

- The `.def` files in the repository serve as a baseline/fallback
- Generated files should match the checked-in versions
- Differences indicate new symbols or build configuration changes
