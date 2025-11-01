# MSVC Library Support for bibfor_ext and asterGC

## Overview

This document describes the implementation of separate `.def` file support for two additional libraries in the MSVC build:
- **bibfor_ext**: Third-party interface functions from `bibfor/third_party_interf/`
- **asterGC**: Garbage collection and memory management from `libs/gc/`

## Background

The MSVC build initially only processed four main libraries:
- `bibfor` - Fortran library (main code with INT64 support)
- `bibc` - C library
- `bibcxx` - C++ library
- `aster` - Main aster library

However, the wscript files actually build two additional separate libraries:
1. **bibfor_ext** (`asterbibfor_ext`) - Built from `bibfor/third_party_interf/*.F90`
   - Contains interfaces to external libraries (MUMPS, PETSc, ELPA)
   - Compiled WITHOUT INT64 to match external library expectations
   
2. **asterGC** (`astergc`) - Built from `libs/gc/*.{F90,cxx}`
   - Contains garbage collection and memory management functions
   - Mixed Fortran/C++ code

## Problem

The original MSVC build system (`msvc/msvc_lib.py`) only generated `.lib` files for the four main libraries. This caused unresolved symbol errors because:

1. `bibfor.def` tried to export symbols from `bibfor_ext` (they don't exist in `bibfor.lib`)
2. `asterGC` library had no `.def` file, so no symbols were exported

## Solution

### 1. Created New DEF Generation Scripts

**`msvc/def_gen_bibfor_ext.py`**
- Scans `build/int64/{debug,release}/bibfor/third_party_interf/*.o` object files
- Extracts Fortran symbols (ending with `_`)
- Generates `msvc/bibfor_ext.def`

**`msvc/def_gen_astergc.py`**
- Scans `build/int64/{debug,release}/libs/gc/*.o` object files  
- Extracts all exported symbols (Fortran and C++)
- Generates `msvc/asterGC.def`

### 2. Updated Build System

**`msvc/msvc_lib.py`**
- Added `asterbibfor_ext` and `astergc` to `LibTask` dataclass
- Updated `extract_main_tasks()` to get task objects for both libraries
- Added def file mappings in `msvclibgen` task
- Updated `run_mvsc_lib_gen()` to:
  - Create `.lib` files for both new libraries
  - Set up proper dependency chains between all libraries
  - Remove `.lib` from shared library task outputs

**`msvc/wscript`**
- Added command-line options for new def files:
  - `--bibfor-ext-def` (default: `msvc/bibfor_ext.def`)
  - `--astergc-def` (default: `msvc/asterGC.def`)
- Updated `configure()` to store paths in environment variables:
  - `self.env.BIBFOR_EXT_DEF`
  - `self.env.ASTERGC_DEF`

**`pixi.toml`**
- Added new tasks under `[feature.dev.target.win-64.tasks]`:
  - `def-bibfor-ext`: Generate bibfor_ext.def
  - `def-astergc`: Generate asterGC.def

### 3. Library Dependencies

The dependency chain is now:
```
aster.lib (asterlib)
  ├─ bibfor.lib (asterbibfor)
  ├─ bibfor_ext.lib (asterbibfor_ext)  [NEW]
  ├─ bibc.lib (asterbibc)
  ├─ bibcxx.lib (asterbibcxx)
  └─ AsterGC.lib (astergc)              [NEW]

bibcxx.lib
  ├─ aster.lib
  ├─ bibc.lib
  ├─ bibfor.lib
  ├─ bibfor_ext.lib                     [NEW]
  └─ AsterGC.lib                        [NEW]

bibc.lib
  ├─ bibcxx.lib
  ├─ bibfor.lib
  └─ bibfor_ext.lib                     [NEW]

bibfor.lib
  ├─ bibcxx.lib
  └─ bibc.lib

bibfor_ext.lib                          [NEW]
  ├─ bibcxx.lib
  └─ bibc.lib

AsterGC.lib                             [NEW]
  ├─ bibc.lib
  └─ bibfor.lib
```

## Usage

### Generating DEF Files

After building the project, regenerate the def files if needed:

```cmd
pixi run def-bibfor-ext
pixi run def-astergc
```

Or regenerate all def files:
```cmd
pixi run def-bibfor
pixi run def-bibfor-ext
pixi run def-bibc
pixi run def-bibcxx
pixi run def-astergc
```

### Build Process

The def files are automatically used during the build:

```cmd
# Debug build
pixi run buildd
pixi run installd

# Release build  
pixi run build
pixi run install
```

### Custom DEF File Paths

You can specify custom paths during configure:

```cmd
python waf configure --bibfor-ext-def=path/to/custom_bibfor_ext.def --astergc-def=path/to/custom_astergc.def
```

## File Structure

```
msvc/
├── def_gen_fc.py              # Generate bibfor.def (existing)
├── def_gen_bibfor_ext.py      # Generate bibfor_ext.def (NEW)
├── def_gen_c.py               # Generate bibc.def (existing)
├── def_gen_cpp.py             # Generate bibcxx.def (existing)
├── def_gen_astergc.py         # Generate asterGC.def (NEW)
├── bibfor.def                 # Main Fortran library exports
├── bibfor_ext.def             # Third-party interface exports (NEW)
├── bibc.def                   # C library exports
├── bibcxx.def                 # C++ library exports
├── asterGC.def                # GC library exports (NEW)
├── aster.def                  # Main aster library exports
├── msvc_lib.py                # MSVC library generation (UPDATED)
└── wscript                    # MSVC build configuration (UPDATED)
```

## Notes

### Why Separate bibfor_ext?

The `bibfor_ext` library is built separately because:
1. It interfaces with external libraries (MUMPS, PETSc, ELPA)
2. External libraries may use different integer sizes (typically int32)
3. Main `bibfor` uses INT64 for better precision
4. Separate compilation allows different `defines` (`WITHOUT_INT64` for bibfor_ext)

### Symbol Counts

Expected symbol counts (approximate):
- **bibfor.def**: ~8000-10000 symbols (main Fortran code)
- **bibfor_ext.def**: ~57 symbols (third-party interfaces)
- **bibc.def**: ~1000-2000 symbols (C code)
- **bibcxx.def**: ~3000-5000 symbols (C++ code)
- **asterGC.def**: ~10-20 symbols (GC functions)
- **aster.def**: ~100-200 symbols (main library)

### Troubleshooting

**Problem**: "No object files found" when running def generation scripts

**Solution**: Build the project first:
```cmd
pixi run buildd  # or pixi run build
```

**Problem**: Unresolved symbols during linking

**Solution**: Regenerate all def files after rebuild:
```cmd
pixi run def-bibfor
pixi run def-bibfor-ext
pixi run def-bibc
pixi run def-bibcxx
pixi run def-astergc
```

**Problem**: Wrong symbols in def file

**Solution**: The def generation scripts use `llvm-nm`. Make sure LLVM tools are in PATH (should be included in pixi environment).

## Implementation Checklist

- [x] Create `msvc/def_gen_bibfor_ext.py`
- [x] Create `msvc/def_gen_astergc.py`
- [x] Update `msvc/msvc_lib.py` - Add bibfor_ext and asterGC to LibTask
- [x] Update `msvc/msvc_lib.py` - Update extract_main_tasks()
- [x] Update `msvc/msvc_lib.py` - Add def file mappings
- [x] Update `msvc/msvc_lib.py` - Update run_mvsc_lib_gen()
- [x] Update `msvc/msvc_lib.py` - Update create_msvclibgen_task()
- [x] Update `msvc/wscript` - Add command-line options
- [x] Update `msvc/wscript` - Store def paths in configure()
- [x] Update `pixi.toml` - Add def-bibfor-ext task
- [x] Update `pixi.toml` - Add def-astergc task
- [x] Create this documentation

## Testing

After implementation, test the build:

1. Clean build directory:
   ```cmd
   rmdir /s /q build
   ```

2. Configure:
   ```cmd
   python waf configure_debug --msvc-entry
   ```

3. Build:
   ```cmd
   pixi run buildd
   ```

4. Verify all libraries exist:
   ```cmd
   dir build\int64\debug\bibfor\bibfor.lib
   dir build\int64\debug\bibfor\bibfor_ext.lib
   dir build\int64\debug\bibc\bibc.lib
   dir build\int64\debug\bibcxx\bibcxx.lib
   dir build\int64\debug\libs\AsterGC.lib
   dir build\int64\debug\bibc\aster.lib
   ```

5. Check for unresolved symbols in build log

## Related Documentation

- `msvc/FORTRAN_UNRESOLVED_REFERENCES.md` - Analysis of the original problem
- `msvc/DEF_AUTO_GENERATION.md` - DEF file generation process
- `msvc/DEF_INTEGRATION_SUMMARY.md` - DEF file integration summary
- `msvc/MSVC_README.md` - General MSVC build documentation

## References

- Code Aster MSVC Guidelines: `.github/copilot-instructions.md`
- WAF Build System: https://waf.io/
- MSVC DEF Files: https://learn.microsoft.com/en-us/cpp/build/reference/module-definition-dot-def-files

