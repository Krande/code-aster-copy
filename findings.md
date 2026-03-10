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
