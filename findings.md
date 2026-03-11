# Root Cause Analysis: code_aster Windows Test Failures

**Date:** 2026-03-07 to 2026-03-09
**Branch:** win-support

## Executive Summary

**ROOT CAUSE: Jeveux COMMON blocks not shared between DLLs on Windows**

Upstream commit `245fa8f519` ("[#33331] do not enable '-i8' option in third_party_interf",
Apr 2025) split `bibfor.dll` into two DLLs:
- `bibfor.dll`: main Fortran code (compiled with `-i8`)
- `bibfor_ext.dll`: third_party_interf/ code including `amumph.F90` (compiled WITHOUT `-i8`)

On Linux, shared libraries share COMMON blocks through the global symbol table, so
Jeveux arrays like `zk24` in the `kvarje` COMMON block are shared across all `.so` files.
**On Windows, each DLL gets its own private copy of COMMON blocks.** This means:

1. `jeveuo` (in `bibfor.dll`) computes index `jrefa` relative to `bibfor.dll`'s copy of `zk24`
2. `amumph` (in `bibfor_ext.dll`) uses `zk24(jrefa-1+11)` with `bibfor_ext.dll`'s copy of `zk24`
3. Since the two `zk24` arrays are at different memory addresses, the index is wrong → **Access Violation**

**Evidence:**
- Build 7 PDB traces show crash at `amumph.F90:207`: `if (zk24(jrefa-1+11) .eq. 'MATR_DISTR')`
- `.def` files for neither `bibfor.dll` nor `bibfor_ext.dll` export COMMON block symbols
- Old 17.2.8 version had ALL Fortran in a single `bibfor.dll` → no COMMON block sharing issue
- Historical data: 84.69% pass rate (17.2.x, single DLL) → 17.85% (17.3.x, split DLLs)

**Fix (Build 10):** Added `!DEC$ ATTRIBUTES DLLEXPORT/DLLIMPORT` directives to `jeveux.h`
for the 6 Jeveux COMMON blocks, conditioned on `ASTER_PLATFORM_MSVC64` and `WITHOUT_INT64`.

**Result:** Pass rate jumped from **19% → 85%** (1924/2276 passed, 352 failed).

---

## Detailed Timeline

### The Regression
- **Feb 2025** (17.2.x): 84.69% pass rate — single `bibfor.dll` with ALL Fortran
- **Apr 2025**: Upstream commit `245fa8f519` splits `bibfor_ext.dll` from `bibfor.dll`
- **Nov 2025** (17.3.x): 17.85% pass rate — split DLLs, COMMON blocks not shared

### Upstream Commit That Caused It
```
245fa8f519 [#33331] do not enable '-i8' option in third_party_interf
Author: Mathieu Courtois <mathieu.courtois@edf.fr>
Date:   Mon Apr 14 22:02:16 2025 +0200
```

This commit:
- Splits `get_srcs("**/*.F90")` into `src_i8` (main) and `src_ext` (third_party_interf)
- Creates separate `bibfor_ext` shared library target
- Applies `-i8`/`-r8` only to `bibfor` via `use=["INT64"]`, not to `bibfor_ext`

---

## Stack Trace (Build 7 with PDB)

```
Exception Code: 0xC0000005 (Access Violation)
  Attempt to read from address 0x0000025B227065D8

[ 5] AMUMPH    amumph.F90:207        ← CRASH
[ 6] TLDLG3    tldlg3.F90:238
[ 7] PRERE1    prere1.F90:167
[ 8] PRERE2    prere2.F90:84
[ 9] MATRIX_FACTOR  matrix_factor.F90:86
[10] LinearSolver::factorize  LinearSolver.cxx:135
```

**amumph.F90:207:**
```fortran
call jeveuo(matas//'.REFA', 'L', jrefa)      ! line 206: jeveuo in bibfor.dll
if (zk24(jrefa-1+11) .eq. 'MATR_DISTR') then ! line 207: zk24 in bibfor_ext.dll → CRASH
```

---

## Why It Only Affects Windows

On Linux/macOS:
- Shared libraries (`.so`) share a single global symbol table
- COMMON blocks like `kvarje` are resolved to the same memory across all libraries
- `zk24` in `libaster.so`, `libbibfor.so`, and `libbibfor_ext.so` all point to the same memory

On Windows:
- DLLs have isolated symbol tables
- COMMON blocks are private to each DLL unless explicitly exported with `DATA` attribute in `.def`
- `zk24` in `bibfor.dll` and `zk24` in `bibfor_ext.dll` are **different arrays at different addresses**
- The current `.def` files only export function symbols, not COMMON block data symbols

---

## COMMON Blocks Affected

From `bibfor/include/jeveux.h`:
```fortran
common / i4vaje / zi4(1)    ! integer*4 array
common / ivarje / zi(1)     ! integer*8 array
common / rvarje / zr(1)     ! real*8 array
common / cvarje / zc(1)     ! complex*8 array
common / lvarje / zl(1)     ! logical array
common / kvarje / zk8(1), zk16(1), zk24(1), zk32(1), zk80(1)  ! character arrays
```

And from `jeveux_private.h`:
```fortran
common / izonje / lk1zon, jk1zon, liszon, jiszon
common / ienvje / lbis, lois, lols, lor8, loc8
common / ilocje / iloc
```

ALL of these need to be shared between `bibfor.dll` and `bibfor_ext.dll`.

---

## Build Results History

| Build | Config | OK | NOOK | CRASH | Key Change |
|-------|--------|-----|------|-------|------------|
| 1-3 | Various MUMPS configs | 431 | 7 | 1845 | Baseline |
| 6 | ifx MUMPS (verified local) | ~431 | ~6 | ~1845 | MUMPS compiler irrelevant |
| 7 | + LINKFLAGS fix (PDB!) | ~431 | ~6 | ~1845 | Got line numbers! |
| 8 | + zk16(25)→zk16(1) | ~same | ~same | ~same | Didn't fix (expected) |
| 9 | + DATA exports in .def (no DLLIMPORT) | - | - | - | Built but didn't fix (need DLLIMPORT) |
| 10 | + DLLEXPORT/DLLIMPORT in jeveux.h | **1924** | **352** | **~0** | **85% pass rate!** |

---

## Previous Hypotheses (DISPROVEN)

1. **Fortran ABI mismatch (ifx vs flang)** — DISPROVEN: Same crashes with ifx MUMPS
2. **jeveux.h zk16(25) COMMON block misalignment** — DISPROVEN: `jxlocs` uses `loc()` at runtime
3. **MUMPS struct access** — DISPROVEN: Build 7 PDB shows crash is at Jeveux zk24 access
4. **Upstream code regression** — CORRECT but root cause is DLL COMMON block isolation

---

## Fixes Applied (Infrastructure)

### Build 7: PDB Debug Symbols
- **Bug**: `build.bat` set `/DEBUG:FULL` in `LDFLAGS` but waf PDB detection checks `LINKFLAGS`
- **Fix**: Set `/DEBUG:FULL` in both `LINKFLAGS` and `LDFLAGS`
- **Result**: 12 PDB files (bibfor.pdb 46MB, bibfor_ext.pdb 1.1MB, bibcxx.pdb 98MB)

### MUMPS ifx Build (builds 3-6)
- ifx with ILP64, PDB files, debug symbols

### Infrastructure
- Fixed stale repodata.json hashes
- win_stacktrace.c: Enhanced symbol search path
- Added `ninja` to test dependencies

---

## MED File I/O Crashes (~64 crashes, 3.6%)

Separate issue, not yet investigated. May also be COMMON-block-related
if MED interface code is in a different DLL.

---

## Fix Applied (Build 10): DLLEXPORT/DLLIMPORT for Jeveux COMMON Blocks

### What Was Done

Three changes were needed:

**1. `bibfor/include/jeveux.h`** — Added Intel Fortran `!DEC$` directives:
```fortran
#if defined(ASTER_PLATFORM_MSVC64) && defined(WITHOUT_INT64)
!DEC$ ATTRIBUTES DLLIMPORT :: /I4VAJE/, /IVARJE/, /RVARJE/
!DEC$ ATTRIBUTES DLLIMPORT :: /CVARJE/, /LVARJE/, /KVARJE/
#elif defined(ASTER_PLATFORM_MSVC64)
!DEC$ ATTRIBUTES DLLEXPORT :: /I4VAJE/, /IVARJE/, /RVARJE/
!DEC$ ATTRIBUTES DLLEXPORT :: /CVARJE/, /LVARJE/, /KVARJE/
#endif
```
- `bibfor.dll` (no `WITHOUT_INT64`) → DLLEXPORT (owns the COMMON blocks)
- `bibfor_ext.dll` (has `WITHOUT_INT64`) → DLLIMPORT (uses bibfor's copy)
- Non-Windows builds → no change

**2. `bibfor/wscript`** — Added `use=["asterbibfor"]` for bibfor_ext on MSVC:
```python
use_ext = use + ["asterbibfor"] if self.env.ASTER_PLATFORM_MSVC64 else use
```
This ensures `bibfor_ext.dll` links against `bibfor.lib` to resolve the imported COMMON blocks.

**3. `msvc/def_gen_fc.py`** — Added DATA exports for COMMON blocks in the .def generator.

### Why `.def` Alone Wasn't Enough (Build 9)

Intel Fortran on Windows requires `!DEC$ ATTRIBUTES DLLIMPORT` in the **source code** for
data symbols. Without it, the compiler generates code that accesses a local copy of the
COMMON block rather than the imported one via the import address table. A `.def` file only
handles linker-level exports but doesn't inform the compiler to generate the correct
indirection for data imports.

### Result

| Build | Pass | Fail | Pass Rate |
|-------|------|------|-----------|
| 8 (before fix) | 431 | 1845 | 19% |
| **10 (DLLEXPORT/DLLIMPORT)** | **1924** | **352** | **85%** |

The remaining 352 failures are genuine test issues (numerical precision, missing features,
Windows path handling, etc.) rather than the systemic COMMON block crash.

---

## Build 11 Regression: flang-built MUMPS (2026-03-11)

### Summary

**Commit `32eadcc8cb`** switched MUMPS from ifx (Intel Fortran) to flang (LLVM Fortran)
to align with conda-forge's mumps-feedstock, reducing the scope of changes needed for
an upstream PR. This caused the pass rate to **drop from 85% back to ~19%**.

**Test results:** 1843 failing tests out of ~2276
- **1728 (94%)** fail with `MUMPS INFO(1)=-3, INFO(2)=0`
  - 1678 DMUMPS, 58 ZMUMPS, 2 SMUMPS
- **~115 (6%)** fail with other errors (Access Violations, MED I/O, etc.)

### ROOT CAUSE: Fortran ABI Incompatibility (ifx vs flang)

**The `DMUMPS_STRUC` Fortran derived type cannot be shared between compilers.**

code_aster (compiled with ifx) creates instances of `DMUMPS_STRUC` and calls
`CALL DMUMPS(dmpsk)`. MUMPS (compiled with flang) receives this struct and reads
its fields. **The derived type layout differs between ifx and flang** because:

1. **Fortran POINTER descriptors are compiler-specific.** `DMUMPS_STRUC` contains
   dozens of `POINTER` array fields (`A`, `IRN`, `JCN`, `RHS`, `IS`, `S`, etc.).
   Each POINTER is internally represented as a descriptor (dope vector):
   - Intel Fortran: uses Intel's proprietary descriptor format
   - LLVM flang: uses the F18/CFI runtime descriptor format
   - These have **different sizes and layouts**

2. **Even with `SEQUENCE`**, which `DMUMPS_STRUC` uses, the Fortran standard does
   NOT mandate a specific representation for POINTER components. `SEQUENCE` only
   guarantees storage order for non-pointer components.

3. **Field offset corruption**: Since the first POINTER field has a different size
   between compilers, ALL subsequent fields (including `JOB`, `N`, `INFO`, `ICNTL`)
   are at wrong offsets. MUMPS reads garbage values, leading to:
   - `INFO(1)=-3` (workspace allocation error with nonsensical parameters)
   - `INFO(2)=0` (corrupted error detail field)

### What Changed (commit `32eadcc8cb`)

**MUMPS build:**
- Fortran compiler: `ifx` → `flang` (via `${{ compiler('fortran') }}`)
- Integer size: ILP64 (`mumps_int_def64_h.in`) → LP64 (`mumps_int_def32_h.in`)
- BLAS: explicit MKL ILP64 → `find_package(BLAS)` (LP64)
- Build type: `RelWithDebInfo` → `Release` (no PDB)
- Removed Intel-specific flags (`/names:lowercase /assume:underscore /Z7 /traceback`)

**code_aster build.bat:**
- Removed `/integer-size:64 /real-size:64` from global FCFLAGS
- These are now correctly handled by waf via `FCFLAGS_INT64` (only applied to bibfor)

### Evidence

The error pattern is uniform across 94% of failures:
```
** ERROR RETURN ** FROM DMUMPS INFO(1)=   -3
** INFO(2)=               0
```

Stack traces consistently show:
```
bibfor_ext.dll  AMUMPM  amumpm.F90:539   ← code_aster wrapper (ifx)
bibfor_ext.dll  AMUMPD  amumpd.F90:246   ← calls DMUMPS(dmpsk)
                                          ← struct layout mismatch here
→ MUMPS (flang) reads corrupted fields → INFO(1)=-3
```

### Why It Worked Before

Build 10 used ifx for **both** MUMPS and code_aster:
- Same compiler → same POINTER descriptor format → same struct layout
- Same ILP64 flags → consistent integer sizes throughout
- Pass rate: 85% (1924/2276)

### Options to Fix

| Option | Effort | Risk | Conda-forge alignment |
|--------|--------|------|----------------------|
| **A: Build MUMPS with ifx** | Low (revert) | Low (proven) | None - separate recipe |
| **B: Build code_aster with flang** | High | High | Full alignment |
| **C: Use MUMPS C interface** | Medium | Medium | Full alignment |
| **D: ifx variant in conda-forge** | Medium | Low | Partial |

**Option A** (revert to ifx MUMPS) is the immediate fix to restore 85% pass rate.

**Option B** (flang for everything) requires:
- Adding flang INT64 flags (`-fdefault-integer-8`, `-fdefault-real-8`) to waf
- Replacing `!DEC$ ATTRIBUTES DLLEXPORT/DLLIMPORT` (Intel-specific) with
  a portable mechanism for COMMON block sharing (`.def` files + linker flags)
- Significant testing

**Option C** (MUMPS C interface) would modify amump*.F90 to call the C wrapper
(`dmumps_c()` which takes `DMUMPS_STRUC_C`, a C struct) instead of the Fortran
interface. The C struct layout is compiler-independent. This requires refactoring
all MUMPS wrapper files but is the most robust long-term solution.

**Option D** (conda-forge ifx variant) would add an `ifx` variant to the
mumps-feedstock, usable when the consumer needs ifx compatibility.

### Recommendation

Short-term: **Option A** (revert to ifx MUMPS) to restore functionality.

Long-term: **Option C** (MUMPS C interface) provides the best cross-compiler
compatibility and would allow using conda-forge's flang-built MUMPS without issues.

---

## Build 12: Name-Mangling Fix for ifx MUMPS (2026-03-11)

### Problem: ifx MUMPS still showed 94% INFO(1)=-3 failures

After rebuilding MUMPS with ifx (Option A), the same `INFO(1)=-3, INFO(2)=0` failures
persisted at ~94%. Root cause analysis by inspecting the DLL's exported symbols revealed:

**The CMake ifx detection didn't fire.** The condition:
```cmake
if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" AND WIN32)
```
did NOT match for the conda-packaged ifx. As a result, `/names:lowercase /assume:underscore`
was never applied to the MUMPS Fortran compilation.

### Root cause: Symbol Name Mismatch → C Wrapper Called Instead of Fortran Entry

Evidence from `dmumps.dll` symbol inspection:
- Fortran symbols were **UPPERCASE** (`DMUMPS`, `DMUMPS_ANA_DRIVER`, etc.)
- Only `dmumps_` (lowercase+underscore) was the **C wrapper** from `mumps_c.c`

The calling chain with the bug:
1. code_aster: `call dmumps(struct)` → ifx generates reference to `dmumps_` (lowercase+underscore)
2. Linker resolves to `dmumps_` in MUMPS DLL → this is the **C wrapper**, not the Fortran entry
3. C wrapper expects `DMUMPS_STRUC_C` (C struct: raw `double*` pointers = 8 bytes each)
4. But receives `DMUMPS_STRUC` (Fortran struct: full POINTER descriptors ≈ 48 bytes each)
5. Field offsets completely wrong → N reads as 0 → `INFO(1)=-3`

### Fix

Bypassed CMake's unreliable compiler ID detection. Fortran flags are now passed directly
from `build-mumps.bat` via `-DCMAKE_Fortran_FLAGS`:

```bat
set "EXTRA_FC_FLAGS="
where ifx >nul 2>nul
if not errorlevel 1 set "EXTRA_FC_FLAGS=/names:lowercase /assume:underscore"
cmake ... -DCMAKE_Fortran_FLAGS="%EXTRA_FC_FLAGS%" ...
```

Also added `-DMUMPS_USE_IFX=ON` to handle CRT library linkage (ifx's Fortran linker
doesn't automatically include C runtime import libraries when linking C objects into
a Fortran DLL).

### Verification

After fix, DLL symbol inspection confirmed:
- `DMUMPS\x00`: 0 occurrences (was 1 before fix)
- `dmumps_\x00`: 1 occurrence (Fortran entry, now lowercase+underscore)
- `dmumps_ana_driver_`: present (lowercase, confirming flags applied)
This is also the approach other projects use when mixing Fortran compilers.

---

## Baseline Verification (2026-03-11)

### Confirmed: commit `00eccd270b` reproduces 84% pass rate

Built in a git worktree at `C:\Work\code\code-aster-src-working` (detached HEAD at `00eccd270b`).

**Result:** 84% tests passed, 353 tests failed out of 2276 (1923 passed).

This confirms the ILP64 approach (ifx MUMPS + global `/integer-size:64`) as the last known working configuration.

### Configuration at that commit

**MUMPS:**
- `mumps_int_def64_h.in` (MUMPS_INT = INTEGER(8))
- `-DUSE_ILP64:BOOL=ON` → `/integer-size:64` in MUMPS Fortran flags + `-DINTSIZE64` C define
- MKL ILP64 libraries (`mkl_intel_ilp64_dll`, `mkl_intel_thread_dll`, `mkl_core_dll`, `libiomp5md`)
- `-DCMAKE_Fortran_COMPILER=ifx` explicit
- `RelWithDebInfo` build type

**code_aster:**
- Global FCFLAGS had `/integer-size:64 /real-size:64` (ALL Fortran code uses 8-byte integers)
- waf probe: `sizeof(dmpsk%n)` = 8 → `ASTER_MUMPS_INT_SIZE=8`
- MKL LP64 (`mkl_intel_lp64_dll`) for code_aster's own BLAS calls

### Next step: LP64 approach (matching Linux)

Linux conda-forge uses LP64 MUMPS (`mumps_int_def32_h.in`, no `-fdefault-integer-8`).
code_aster's Linux feedstock does NOT set global integer size flags — waf handles via
`FCFLAGS_INT64` (applied only to bibfor, not bibfor_ext). `ASTER_MUMPS_INT_SIZE=4` on Linux.

The current branch is being updated to match this LP64 approach:
- MUMPS: LP64, `mumps_int_def32_h.in`, no `/integer-size:64`
- code_aster: no global `/integer-size:64`, waf-managed FCFLAGS_INT64
- Keep: name-mangling fix, CRT lib fix, debug symbols (`RelWithDebInfo` + `/DEBUG:FULL`)
- Keep: dual variant support (ifx+MKL / flang+netlib)

If the LP64 approach achieves the same ~84% pass rate, it would be the preferred
configuration as it aligns with Linux and avoids global integer size pollution.
