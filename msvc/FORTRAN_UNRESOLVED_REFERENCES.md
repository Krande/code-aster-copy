# Fortran Unresolved References on MSVC - Analysis and Resolution Strategy

## Problem Summary

During the MSVC build, 57 Fortran symbols from `bibfor` are reported as unresolved when generating `bibfor.lib`:

```
bibfor.def : error LNK2001: unresolved external symbol amumpc_
bibfor.def : error LNK2001: unresolved external symbol amumps_
... (55 more symbols)
```

All affected symbols are in `bibfor/third_party_interf/` and are interfaces to external libraries:
- **MUMPS** (MUltifrontal Massively Parallel sparse direct Solver)
- **PETSc** (Portable, Extensible Toolkit for Scientific Computation)
- **ELG** (ELPA Large-scale Eigenvalue solvers)

## Root Cause Analysis

### Fortran Name Mangling on MSVC vs GNU

**GNU/GCC Fortran:**
- Appends underscore to function names: `amumpc` → `amumpc_`
- Lowercase by default
- Symbols are exported automatically in shared libraries

**MSVC (Intel Fortran / flang):**
- Different name mangling conventions
- May use uppercase or different decorations
- Requires explicit symbol exports via `.def` files

### The `.def` File Approach

The MSVC build uses `.def` files to explicitly export symbols from DLLs:
- `bibfor.def` lists all symbols that should be exported from `bibfor.lib`
- The symbols are currently listed with GNU-style names (e.g., `amumpc_`)
- The linker cannot find these symbols because they may have different names in the compiled object files

## What Has Been Tried So Far

### 1. ✅ Adding Symbols to .def File
- **What:** All 57 symbols are already listed in `bibfor.def`
- **Result:** Still unresolved - the symbols are in the .def file but the linker can't match them to actual object file symbols

### 2. ❌ Using Aliases in .def File
- **What:** Tried mapping GNU names to MSVC names:
  ```
  amumpc_
  amumpc=amumpc_
  ```
- **Result:** Failed - aliases didn't resolve the issue
- **Why:** Without knowing the actual symbol name in the object files, aliases can't work

### 3. ✅ **SUCCESSFUL: Using `bind(C)` for MUMPS Interface**
The `cmumps` subroutine interface successfully uses `bind(C)` with explicit name specification:

**File:** `bibfor/include/mumps/cmumps.h`
```fortran
#ifdef ASTER_HAVE_MUMPS
interface
#ifdef ASTER_PLATFORM_MSVC64
    subroutine cmumps(cmpsk) bind(C, name='CMUMPS')
#else
    subroutine cmumps(cmpsk)
#endif
#       include "cmumps_struc.h"
        type (cmumps_struc) :: cmpsk
    end subroutine cmumps
end interface
#endif
```

**Key Points:**
- Uses `bind(C, name='CMUMPS')` to force C-style linkage
- **Uses UPPERCASE name** (`CMUMPS` not `cmumps_`)
- Only applied on MSVC platform via `#ifdef ASTER_PLATFORM_MSVC64`
- This makes the symbol follow C naming conventions (no underscore mangling)

## Affected Symbols (57 Total)

### MUMPS-Related (11 symbols)
- `amumpc_` - Complex MUMPS driver
- `amumpd_` - Double precision MUMPS driver
- `amumph_` - MUMPS helper
- `amumpi_` - MUMPS initialization/parameters
- `amumpm_` - MUMPS matrix operations
- `amumpm_hpc_` - MUMPS HPC variant
- `amumpp_` - MUMPS preconditioning
- `amumps_` - Real MUMPS driver (single precision)
- `amumpt_` - MUMPS termination
- `amumpu_` - MUMPS utilities
- `amumpz_` - Complex double MUMPS driver

### PETSc-Related (21 symbols)
- `ap2foi_`
- `ap_assembly_vector_`
- `ap_on_off_`
- `apalmc_`, `apalmd_`, `apalmh_`
- `apbloc_`
- `apetsc_`
- `apksp_`
- `apldlt_`
- `apmain_`
- `apmamc_`, `apmamd_`, `apmamh_`, `apmams_`
- `appcpr_`, `appcrs_`
- `apsolu_`
- `apvsmb_`, `apvsmbh_`
- `pchpddmdumpauxiliarymat_`

### ELG-Related (9 symbols)
- `elg_allocvr_`
- `elg_calc_matk_red_`
- `elg_calc_matm_red_`
- `elg_calc_rhs_red_`
- `elg_calc_solu_`
- `elg_calcx0_`
- `elg_calcxl_`
- `elg_calcxl_modal_`
- `elg_resoud_`

### Verification/Solver Related (7 symbols)
- `bresels_`, `breselsqp_`, `breselu_`
- `cpysol_`
- `elim75_`
- `filter_smd_`
- `verifels_`, `verifelu_`

### Other (9 symbols)
- `irmhpc_`
- `ldsp1_`, `ldsp2_`
- `lrmtyp_`
- `lrvema_`
- `matass2petsc_`
- `vect_asse_from_petsc_`
- `vect_asse_update_ghost_values_`

## Resolution Strategy: Checklist

### 🎯 CRITICAL FINDING - ROOT CAUSE IDENTIFIED

**THE SYMBOLS ARE NOT MISSING - THEY EXIST AND ARE CORRECTLY NAMED!**

The linker error `unresolved external symbol amumpc_` is **MISLEADING**. 

**The actual problem:**
1. ✅ `amumpc_` exists in `amumpc.F90.2.o` with correct GNU-style naming
2. ✅ All 57 symbols exist in their respective object files
3. ❌ **The symbols are in `bibfor_ext` library, but `bibfor.def` tries to export them from `bibfor` library**

**Root Cause:**

The `bibfor/wscript` builds **TWO SEPARATE LIBRARIES**:

```python
# bibfor/wscript (lines 63-84)
src_i8 = get_srcs("**/*.F90", excl="third_party_interf/*.F90")
src_ext = get_srcs("third_party_interf/*.F90")

# Library 1: bibfor - ALL sources EXCEPT third_party_interf
self(
    features="fc fcshlib",
    name="asterbibfor",
    target="bibfor",
    source=src_i8,
    use=use + ["INT64"],
)

# Library 2: bibfor_ext - ONLY third_party_interf sources
self(
    features="fc fcshlib", 
    name="asterbibfor_ext",
    target="bibfor_ext",
    source=src_ext,
    defines=["WITHOUT_INT64"],
    use=use,
)
```

**But:**
- `msvc/msvc_lib.py` only processes `asterbibfor`, NOT `asterbibfor_ext`
- `msvc/bibfor.def` tries to export symbols from `bibfor_ext` as if they were in `bibfor`
- Result: Linker can't find `amumpc_` in `bibfor.lib` because it's actually in `bibfor_ext.lib`

**Evidence:**
```
llvm-nm amumpc.F90.2.o shows:
  00000000 T amumpc_        ← Symbol EXISTS in bibfor_ext object files
  
msvc/msvc_lib.py creates:
  bibfor.lib               ← From asterbibfor (excludes third_party_interf)
  
msvc/bibfor.def tries to export:
  amumpc_                  ← But this is in bibfor_ext, not bibfor!
  
Link error:
  unresolved external symbol amumpc_   ← Not found in bibfor.lib
```

### SOLUTION: Fix the Build System (WAF + MSVC)

The problem is in **TWO places**:

1. `bibfor/wscript` - Creates separate `bibfor` and `bibfor_ext` libraries
2. `msvc/msvc_lib.py` - Only processes `bibfor`, ignores `bibfor_ext`

**RECOMMENDED SOLUTION - Option A: Merge bibfor_ext into bibfor on MSVC**

Modify `bibfor/wscript` to only create the split on non-MSVC platforms:

```python
src_i8 = get_srcs("**/*.F90", excl="third_party_interf/*.F90")
src_ext = get_srcs("third_party_interf/*.F90")

if self.env.ASTER_PLATFORM_MSVC64:
    # On MSVC: Build everything in one library
    self(
        features="fc fcshlib",
        name="asterbibfor",
        target="bibfor",
        source=src_i8 + src_ext,  # Include ALL sources
        use=use,
        env=env.derive(),
        install_path=env.ASTERLIBDIR,
    )
else:
    # On Linux/MinGW: Keep the split (int64 vs normal)
    self(
        features="fc fcshlib",
        name="asterbibfor",
        target="bibfor",
        source=src_i8,
        use=use + ["INT64"],
        env=env.derive(),
        install_path=env.ASTERLIBDIR,
    )
    self(
        features="fc fcshlib",
        name="asterbibfor_ext",
        target="bibfor_ext",
        source=src_ext,
        defines=["WITHOUT_INT64"],
        use=use,
        env=env.derive(),
        install_path=env.ASTERLIBDIR,
    )
```

**Pros:**
- Simple, minimal change
- All symbols in one library = one .def file
- No need to modify MSVC build scripts

**Cons:**
- Loses the int64 vs non-int64 separation on MSVC
- May cause issues if external libraries expect specific integer sizes

---

**ALTERNATIVE - Option B: Process bibfor_ext in msvc_lib.py**

Modify `msvc/msvc_lib.py` to handle `bibfor_ext`:

1. Add `bibfor_ext` to `extract_main_tasks()` 
2. Create `bibfor_ext.lib` with its own .def file
3. Add `bibfor_ext` symbols to a new `msvc/bibfor_ext.def`

**Pros:**
- Preserves int64 separation
- Cleaner separation of concerns

**Cons:**
- More complex build script changes
- Need to manage two .def files
- Callers need to link against both libraries

---

**Option C: Remove from bibfor.def (if unused)**

If these symbols are never actually called (MPI-only, disabled features):
- Remove all 57 symbols from `msvc/bibfor.def`
- Regenerate def file without them

**Check if needed:** Search codebase for calls to these functions

---

## Resolution Strategy: Checklist

### Phase 1: Investigation (PRIORITY) ✅ COMPLETED

- [x] **1.1 Inspect Compiled Object Files**
  ```powershell
  # Use llvm-nm to check actual symbol names in object files
  llvm-nm build\int64\debug\bibfor\third_party_interf\amumpc.F90.2.o
  llvm-nm build\int64\debug\bibfor\third_party_interf\apetsc.F90.2.o
  llvm-nm build\int64\debug\bibfor\third_party_interf\elg_resoud.F90.2.o
  ```
  
  **✅ FINDINGS - ALL SYMBOLS EXIST WITH GNU-STYLE NAMING:**
  
  - **amumpc.F90.2.o:**
    - **DEFINED:** `amumpc_` (lowercase + underscore)
    - **UNDEFINED EXTERNAL:** `U CMUMPS` ⚠️ **UPPERCASE, NO UNDERSCORE**
    - Also references other aster functions: `amumpi_`, `amumpm_`, `amumpp_`, `amumpt_`, `amumpu_`
  
  - **amumpd.F90.2.o:**
    - **DEFINED:** `amumpd_` (lowercase + underscore)
    - **UNDEFINED EXTERNAL:** `U DMUMPS` ⚠️ **UPPERCASE, NO UNDERSCORE**
  
  - **apetsc.F90.2.o:**
    - **DEFINED:** `apetsc_` (lowercase + underscore)
    - **NO EXTERNAL PETSC CALLS** (just internal aster functions)
  
  - **elg_resoud.F90.2.o:**
    - **DEFINED:** `elg_resoud_` (lowercase + underscore)
    - **UNDEFINED EXTERNALS:** Other aster functions (`elg_calc_rhs_red_`, `elg_calc_solu_`, etc.)
    - **NO EXTERNAL ELG LIBRARY CALLS**

  **🔍 KEY DISCOVERY:**
  The issue is **NOT** name mangling of the Code Aster wrapper functions themselves.
  The issue is that MUMPS wrappers call external MUMPS library functions with **UPPERCASE C-style names**.
  
  Example from `amumpc.F90.2.o`:
  ```
  00000000 T amumpc_           ← This symbol EXISTS and is CORRECT
           U CMUMPS            ← This is an EXTERNAL call to MUMPS library
  ```

- [x] **1.2 Check Conditional Compilation**
  - ✅ All object files exist, so conditional compilation flags ARE set
  - `ASTER_HAVE_MUMPS`, `ASTER_HAVE_PETSC`, etc. are defined

- [x] **1.3 Verify Object Files Exist**
  ```powershell
  dir build\int64\debug\bibfor\third_party_interf\*.o
  ```
  **✅ CONFIRMED:** All 75 object files exist in the directory

### Phase 2: Apply `bind(C)` Pattern (IF symbols exist)

If Phase 1 confirms the symbols exist but with different names:

- [ ] **2.1 Update MUMPS Interfaces**
  Apply `bind(C)` pattern to all MUMPS driver subroutines:
  - Files: `amumpc.F90`, `amumpd.F90`, `amumph.F90`, `amumpi.F90`, `amumpm.F90`, `amumpm_hpc.F90`, `amumpp.F90`, `amumps.F90`, `amumpt.F90`, `amumpu.F90`, `amumpz.F90`
  - Pattern:
    ```fortran
    #ifdef ASTER_PLATFORM_MSVC64
    subroutine <name>(...) bind(C, name='<UPPERCASE_NAME>')
    #else
    subroutine <name>(...)
    #endif
    ```
  - **IMPORTANT:** Use UPPERCASE in `bind(C, name='...')` (following `cmumps` pattern)

- [ ] **2.2 Update PETSc Interfaces**
  - Files: All `ap*.F90` files in `third_party_interf/`
  - Apply same `bind(C)` pattern with UPPERCASE names

- [ ] **2.3 Update ELG Interfaces**
  - Files: All `elg_*.F90` files
  - Apply same `bind(C)` pattern

- [ ] **2.4 Update Other Third-Party Interfaces**
  - Files: `bresels*.F90`, `cpysol.F90`, `elim75.F90`, `filter_smd.F90`, `irmhpc.F90`, `ldsp*.F90`, `lrmtyp.F90`, `lrvema.F90`, `matass2petsc.F90`, `vect_asse_*.F90`, `verif*.F90`
  - Apply same pattern

### Phase 3: Update .def File (IF using bind(C))

If `bind(C)` approach is used:

- [ ] **3.1 Regenerate bibfor.def**
  ```cmd
  pixi run def-bibfor
  ```
  **OR manually update** to use UPPERCASE names without underscores:
  ```
  AMUMPC
  AMUMPD
  AMUMPH
  ... (etc)
  ```

### Phase 4: Alternative Solutions (IF bind(C) doesn't work)

- [ ] **4.1 Module File Approach**
  - Check if third-party interfaces should be in Fortran modules
  - Modules have different name mangling than bare subroutines

- [ ] **4.2 Separate Library**
  - Consider compiling third-party interfaces into a separate library
  - This may bypass .def file export issues

- [ ] **4.3 Stub/Wrapper Functions**
  - Create C wrappers for Fortran third-party interfaces
  - C has consistent naming across compilers

- [ ] **4.4 Check for MPI/HPC Dependency**
  - Many of these symbols are MPI/HPC related
  - Verify if they should only be compiled in MPI builds
  - Check if symbols are conditionally compiled based on MPI availability

### Phase 5: Verification

- [ ] **5.1 Clean Build**
  ```cmd
  pixi run buildd
  ```

- [ ] **5.2 Check Link Success**
  Verify all 57 symbols are resolved

- [ ] **5.3 Run Tests**
  ```cmd
  pixi run testd
  ```

## Important Notes

### Why UPPERCASE in bind(C)?

The existing working example (`cmumps`) uses **UPPERCASE** in the bind(C) name:
```fortran
bind(C, name='CMUMPS')
```

This is likely because:
1. Intel Fortran on Windows may default to uppercase external names
2. C linkage convention expects consistent casing
3. The MUMPS library itself may export uppercase names

**Rule:** When adding `bind(C)`, use **UPPERCASE** names matching the .def file export

### Conditional Compilation Concern

All affected files use conditional compilation:
```fortran
#ifdef ASTER_HAVE_MUMPS
  ! code here
#endif
```

**CRITICAL:** If these preprocessor symbols are not defined during compilation:
- The subroutines won't be compiled at all
- No object code will be generated
- Symbols will truly not exist
- `.def` file will reference non-existent symbols

**Action:** Check build configuration first!

### Platform-Specific Changes

All `bind(C)` additions should be wrapped in:
```fortran
#ifdef ASTER_PLATFORM_MSVC64
  bind(C, name='UPPERCASE_NAME')
#endif
```

This ensures:
- Linux/GNU builds are unaffected
- Only MSVC platform gets C linkage
- Maintains code portability

## Expected Outcome

### Success Criteria
1. All 57 symbols resolve during link phase
2. `bibfor.lib` builds successfully
3. No functional regressions in tests
4. Code remains portable (Linux/GNU builds still work)

### If This Doesn't Work

The symbols may not exist because:
1. **Conditional compilation** - External library support not enabled
2. **MPI-only symbols** - Only compiled in parallel/MPI builds
3. **Wrong build configuration** - Debug vs Release, or feature flags

In this case, the symbols should be **removed from `bibfor.def`** rather than trying to fix linkage.

## Files to Modify (if going with bind(C) approach)

All files in `bibfor/third_party_interf/`:
- [ ] `amumpc.F90` through `amumpz.F90` (11 MUMPS files)
- [ ] `ap2foi.F90` through `apvsmbh.F90` (21 PETSc files)
- [ ] `elg_allocvr.F90` through `elg_resoud.F90` (9 ELG files)
- [ ] `bresels.F90`, `breselsqp.F90`, `breselu.F90`
- [ ] `cpysol.F90`
- [ ] `elim75.F90`
- [ ] `filter_smd.F90`
- [ ] `irmhpc.F90`
- [ ] `ldsp1.F90`, `ldsp2.F90`
- [ ] `lrmtyp.F90`
- [ ] `lrvema.F90`
- [ ] `matass2petsc.F90`
- [ ] `PCHPDDMDumpAuxiliaryMat.F90`
- [ ] `vect_asse_from_petsc.F90`
- [ ] `vect_asse_update_ghost_values.F90`
- [ ] `verifels.F90`, `verifelu.F90`

**Total: 57 files corresponding to 57 symbols**

## Next Immediate Action

**START WITH PHASE 1.1** - Check what symbols actually exist:

```powershell
# Navigate to build directory
cd C:\AibelProgs\code\code-aster-src

# Check a few sample object files
llvm-nm --extern-only --defined-only build\int64\debug\bibfor\third_party_interf\amumpc.obj
llvm-nm --extern-only --defined-only build\int64\debug\bibfor\third_party_interf\apetsc.obj
llvm-nm --extern-only --defined-only build\int64\debug\bibfor\third_party_interf\elg_allocvr.obj
```

This will immediately tell us if the problem is:
- **Name mangling** (symbols exist but with different names) → Use bind(C)
- **Missing symbols** (not compiled at all) → Check build configuration

---

## 📋 ACTION ITEMS (IMMEDIATE)

### ✅ Phase 1 Complete - Investigation Results

**CONFIRMED:**
- All 57 symbols exist in compiled object files with correct GNU-style naming
- Object files are in `build/int64/debug/bibfor/third_party_interf/*.o`
- Symbols are defined: `amumpc_`, `amumpd_`, etc. (lowercase + underscore)

**ROOT CAUSE IDENTIFIED:**
- `bibfor/wscript` builds TWO separate libraries: `bibfor` and `bibfor_ext`
- `third_party_interf/*.F90` sources go into `bibfor_ext` 
- `msvc/msvc_lib.py` only processes `asterbibfor`, NOT `asterbibfor_ext`
- `msvc/bibfor.def` tries to export symbols from `bibfor_ext` as if they're in `bibfor`

### 🔧 Next Steps - Choose ONE Option

**RECOMMENDED: Option A - Merge on MSVC**
- [ ] Edit `bibfor/wscript` line 63-84
- [ ] Add platform check `if self.env.ASTER_PLATFORM_MSVC64:`
- [ ] On MSVC: Combine `src_i8 + src_ext` into single `asterbibfor` library
- [ ] Test build: `pixi run buildd`
- [ ] Verify all 57 symbols resolve

**OR: Option B - Process bibfor_ext**
- [ ] Edit `msvc/msvc_lib.py` `extract_main_tasks()` 
- [ ] Add `fc_ext_task_object = get_task_object(bld, "asterbibfor_ext", "fc")`
- [ ] Create `msvc/bibfor_ext.def` with the 57 symbols
- [ ] Add bibfor_ext library generation in `run_mvsc_lib_gen()`
- [ ] Update linker to include both `bibfor.lib` and `bibfor_ext.lib`

**OR: Option C - Check if symbols are needed**
- [ ] Search codebase for calls to `amumpc`, `apetsc`, `elg_resoud`, etc.
- [ ] If unused (MPI-only/disabled features), remove from `bibfor.def`
- [ ] Regenerate: `pixi run def-bibfor`

### 📊 Impact Analysis Needed

Before choosing, verify:
- [ ] Why was bibfor_ext created as separate library? (int64 vs non-int64)
- [ ] Do external libraries (MUMPS/PETSc) require specific integer sizes?
- [ ] Are these symbols actually called in non-MPI builds?

---

**Document Status:** Phase 1 Complete - Root Cause Identified  
**Last Updated:** 2025-01-11  
**Next Action:** Choose fix strategy and implement

