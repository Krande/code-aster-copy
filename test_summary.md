# code_aster Windows Test Summary

**Date:** 2026-03-07
**Branch:** win-support
**Package:** code-aster-17.3.12.1-py312_nompi_release_INT64_h94befa5_1
**Platform:** Windows 11 Pro (win-64), Intel Fortran (ifx), MSVC (cl), Python 3.12

## Overall Results

| Metric | Count |
|--------|-------|
| Total tests | 2276 |
| Passed (ctest) | 431 (19%) |
| Failed (ctest) | 1845 (81%) |
| Total runtime | 1111 sec (~18.5 min) |

### Breakdown by Diagnostic (from .mess files)

| Diagnostic | Count | Description |
|------------|-------|-------------|
| `<F>_ABNORMAL_ABORT` | 1820 | Fatal crash (mostly Access Violation) |
| `NO_TEST_RESU` | 13 | Ran without crash but no validation |
| `NOOK_TEST_RESU` | 7 | Ran but validation failed |
| `<A>_ALARM` | 3 | Completed with warnings |
| `OK` | 2 | Fully passed |
| `<S>_ERROR` | 2 | Expected error tests |

> Note: 431 tests passed at the ctest level, but only ~1845 produced `.mess` files. The
> remaining ~431 passed tests either don't write to the nompi output directory or completed
> without triggering the crashing code paths.

---

## Root Cause #1: MUMPS Access Violation (95.6% of crashes)

**1716 out of 1795 Access Violation crashes** occur in `amumph_` during matrix factorization.

### Stack Trace

```
[ 5] amumph_                  (bibfor_ext)   <-- MUMPS Fortran interface
[ 6] tldlg3_                  (bibfor)
[ 7] prere1_                  (bibfor)
[ 8] prere2_                  (bibfor)
[ 9] matrix_factor_           (bibfor)       <-- Matrix factorization entry
[10] preres_                  (bibfor)
[11] nmprac_                  (bibfor)
[12] accel0_                  (bibfor)
[13] nminit_                  (bibfor)
[14] op0070_                  (bibfor)       <-- STAT_NON_LINE / DYNA_NON_LINE
```

**Exception:** `Attempt to read from address 0x000002930C5007D8`
**Crash address:** `0x00007FFEC98BCC18` (consistent across all 1716 crashes)

### Hypothesis

The crash is in the MUMPS solver interface (`amumph_`), which wraps calls to the MUMPS
library (`mumps-seq` conda package). Likely causes:
- **Integer size mismatch**: code_aster is built with INT64 (`-DASTER_INT_SIZE=8`) but
  MUMPS may expect 32-bit integers for certain struct fields
- **Calling convention mismatch**: The MUMPS Fortran `bind(C)` removal in the previous
  session may have introduced issues, or the mumps-seq package itself has symbol naming
  problems
- **MUMPS struct layout**: The `zmumps_struc.h` / `smumps_struc.h` Fortran derived types
  must match the compiled MUMPS library's struct layout exactly

---

## Root Cause #2: MED File I/O Crashes (3.6% of crashes)

**~64 crashes** in MED library operations, split into three sub-patterns:

### 2a. MED File Open via medfwrap (25 crashes)
```
mfivop_ (medfwrap) -> as_mfivop_ -> ... -> op0039_ (IMPR_RESU)
```
**Exception:** `Attempt to write to address 0x00000000000000C8` (near-NULL pointer)

### 2b. MED C String Handling (25 crashes)
```
MED2cstring (medC) -> MFIFVOP -> mfivop_ (medfwrap) -> ...
```
**Exception:** `Attempt to read from address 0x00007FFEC0ED244B`

### 2c. MED Type Reading (14 crashes)
```
lrmtyp_ (bibfor_ext) -> get_med_types_ -> eval_formula_ (bibcxx)
```

### Hypothesis

String passing between Fortran and C in MED library. The `fix_exports.py` symbol aliasing
ensures linkage works, but the actual **string descriptor** format may differ. Intel Fortran
passes strings with hidden length arguments; the C side (`medC`) may not be receiving them
correctly. The near-NULL pointer write (`0x000000C8`) suggests a struct field offset is
being used as a pointer.

---

## Root Cause #3: Jeveux/MFront Memory Crashes (4 crashes)

```
jedema_ (bibfor) -> etenca_ -> debcal_ -> calcul_ -> mgis_debug_ (bibcxx)
```
**Exception:** `Attempt to read from address 0xFFFFFFFFFFFFFFFF` (uninitialized pointer, -1)

Likely an MFront/MGIS integration issue where the behavior compilation fails silently and
leaves uninitialized data.

---

## Root Cause #4: MFront Runtime Compilation (popup issue)

**Symptom:** Windows popup "mfront: callCmake: can't build target 'all'"

The `comp012*` tests (6 tests) and other MFront-dependent tests fail because mfront tries
to compile material behavior laws at runtime using CMake + Ninja. The test environment has
`cmake` but **not `ninja`**, so the build fails.

**Fix:** Add `ninja` to test run requirements in `recipe.yaml`.

---

## Root Cause #5: NOOK_TEST_RESU Failures (7 tests)

| Test | Issue |
|------|-------|
| asrun01b | Subprocess returns 127 (command not found) |
| mesh001m | Mesh element count mismatch (28!=35, 7!=22) |
| sdlv126e | Numerical tolerance failure |
| sdlv126f | Numerical tolerance failure |
| ssls134a | Timing values all zero |
| vocab01b | Completely wrong result (9009 vs 0) |
| zzzz130b | Timing values wrong |

---

## Priority for Debugging

1. **MUMPS crash** -- Fixing this single issue would resolve ~1716 of 1845 failures (93%)
2. **MED I/O crashes** -- Would resolve ~64 additional failures
3. **MFront runtime** -- Add ninja to test deps (trivial fix)
4. **NOOK failures** -- Lower priority, likely platform-specific numerical differences

---

## Passing Test Categories

Tests that pass tend to be:
- Simple linear algebra (sdll*, ssll*, sslp*)
- Thermal problems (tpl*, ttl*)
- THM coupled problems (wtn*)
- Utility/infrastructure tests (zzzz*, supv*, func*, table*)
- Tests using MULT_FRONT solver instead of MUMPS
- Tests that don't read/write MED files via the crashing path
