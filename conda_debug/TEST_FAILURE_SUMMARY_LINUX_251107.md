# Code_Aster Linux Conda Build - Test Failure Analysis

**Date:** 2025-11-07  
**Test Results:** 96.18% passing (2188/2275 tests)  
**Failed Tests:** 87 tests (3.82%)

---

## Executive Summary

The Code_Aster conda package for Linux shows excellent test pass rate of **96.18%**. The 87 failing tests have been categorized into specific root causes, with the majority being:
1. **Missing runtime dependencies** (60 tests, 68.97%)
2. **MFront/MGIS configuration issues** (21 tests, 24.14%)
3. **Segmentation faults** (6 tests, 6.90%)

---

## Failure Categories

### 1. Missing Runtime Dependencies (60 tests - 68.97%)

#### 1.1 Missing medcoupling (27 tests - 31.03%)
**Impact:** Highest number of failures  
**Root Cause:** The `medcoupling` Python module is missing from runtime environment. This was in `host` dependencies but not in `run` dependencies in the conda recipe.

**Affected Tests:**
```
zzzz513a, zzzz142a, hplv111c, hplv106a, zzzz503h, zzzz154a, hplv106c, hplv106b, 
zzzz513b, hplv108e, mesh001e, hplv108c, hplv108d, hplv108b, zzzz513e, zzzz143a, 
hplv111a, zzzz509q, hplv111d, zzzz509r, mesh001a, zzzz507a, zzzz513c, zzzz505i, 
hplv111b, zzzz505g, hplv108a
```

**Error Pattern:**
```python
ModuleNotFoundError: No module named 'medcoupling'
```

**Fix:** ✅ Already fixed - Added `medcoupling` to run dependencies

---

#### 1.2 Missing homard (17 tests - 19.54%)
**Root Cause:** HOMARD mesh adaptation tool not available in runtime

**Affected Tests:**
```
zzzz175b, wtnl100f, forma01b, forma01c, zzzz356a, tpll01h, tplp107b, zzzz319b, 
ssnv173k, zzzz319a, zzzz121c, zzzz121d, tplp305d, wtnl100c, forma11a, tpll01j, 
zzzz121e
```

**Error Pattern:**
```
Le fichier homard est inconnu
```

**Fix Required:** Add `homard` to conda dependencies or mark these tests as expected failures

---

#### 1.3 Missing xmgrace (5 tests - 5.75%)
**Root Cause:** XMGrace plotting tool not available

**Affected Tests:**
```
sdnl105b, forma10a, forma10b, forma30b, sdnl105a
```

**Error Pattern:**
```
Le fichier xmgrace n'existe pas.
```

**Fix Required:** Add `xmgrace` to conda dependencies or skip visualization tests

---

#### 1.4 Missing PETSc (3 tests - 3.45%)
**Root Cause:** PETSc solver not compiled/installed in this build

**Affected Tests:**
```
zzzz512g, zzzz512f, zzzz512h
```

**Error Pattern:**
```
<FERMETUR_10> - Le solveur "PETSc" n'est pas installé dans cette version de Code_Aster.
```

**Fix Required:** Either enable PETSc in the build or mark these tests as expected failures

---

#### 1.5 Missing asrun (2 tests - 2.30%)
**Root Cause:** Legacy asrun module not available

**Affected Tests:**
```
ssna110a, ssnv228b
```

**Error Pattern:**
```
ModuleNotFoundError: No module named 'asrun'
```

**Fix Required:** Add asrun module or update tests to not require it

---

#### 1.6 Missing aspell (1 test - 1.15%)
**Root Cause:** Aspell spell checker not available with correct permissions

**Affected Tests:**
```
supv002a
```

**Error Pattern:**
```
PermissionError: [Errno 13] Permission denied: 'aspell'
```

**Fix Required:** Add `aspell` to conda dependencies

---

#### 1.7 Missing MISS3D (1 test - 1.15%)
**Root Cause:** MISS3D soil-structure interaction tool has permission issues

**Affected Tests:**
```
sdls118a
```

**Error Pattern:**
```
sh: 1: run_miss3d: Permission denied
```

**Fix Required:** Fix MISS3D installation and permissions or mark test as expected failure

---

#### 1.8 Missing Fortran Compiler (1 test - 1.15%)
**Root Cause:** Fortran compiler path from build environment is hardcoded and not available at runtime

**Affected Tests:**
```
zzzz409a
```

**Error Pattern:**
```
FileNotFoundError: Fortran compiler not found: /home/krande/.../x86_64-conda-linux-gnu-gfortran
```

**Fix Required:** Update code to use runtime compiler or provide compiler in runtime environment

---

### 2. MFront/MGIS Configuration Issues (21 tests - 24.14%)

#### 2.1 MGIS Configuration Error (20 tests - 22.99%)
**Root Cause:** The `ASTER_BEHAVIOUR_LIB` configuration variable is `None`, causing string concatenation errors when trying to load MFront behavior libraries.

**Affected Tests:**
```
ssnv518a, hsna106a, comp012f, mfron02j, ssnv163e, hsnv135a, ssnv205c, comp012g, 
ssnv232c, ssnv205b, hsnv136d, comp012h, ssnv520a, comp012i, zzzz509f, wtnv122e, 
ssnv519a, ssnv163d, ssnv181b, wtnv122d
```

**Error Pattern:**
```python
TypeError: can only concatenate str (not "NoneType") to str
File ".../MGISHelper.py", line 47, in from_embedded
    os.environ["ASTER_LIBDIR"], "lib" + config["ASTER_BEHAVIOUR_LIB"] + ".so"
```

**Fix Required:** Set `ASTER_BEHAVIOUR_LIB` properly in configuration during build/install

---

#### 2.2 MFront Compilation Failed (1 test - 1.15%)
**Root Cause:** MFront cannot compile behavior library at runtime

**Affected Tests:**
```
ssnv187x
```

**Error Pattern:**
```
<MFRONT_4> - Le fichier de sortie de MFront libBehaviour.so n'a pas pu être produit.
```

**Fix Required:** Ensure MFront toolchain is properly configured for runtime compilation

---

### 3. Segmentation Faults (6 tests - 6.90%)
**Root Cause:** Memory corruption or library incompatibility issues causing crashes

**Affected Tests:**
```
sslx101a, zzzz413a, ssnl123a, zzzz353a, forma40b, ssnl504a
```

**Error Pattern:**
```
Segmentation fault (core dumped)
```

**Fix Required:** Requires debugging to identify:
- Memory corruption issues
- Library version incompatibilities
- Fortran/C++ ABI issues
- Stack overflow problems

---

### 4. Test Validation Failures (2 tests - 2.30%)

#### 4.1 NOOK_TEST_RESU (2 tests - 2.30%)
**Root Cause:** Test results don't match expected values

**Affected Tests:**
```
ssls150a, sdll123d
```

**Error Pattern:**
```
DIAGNOSTIC JOB : NOOK_TEST_RESU
```

**Fix Required:** Investigate numerical differences or update reference values

---

### 5. Numerical/Solver Issues (1 test - 1.15%)

#### 5.1 Singular Matrix (1 test - 1.15%)
**Root Cause:** Matrix singularity in mechanical calculation

**Affected Tests:**
```
sdls505a
```

**Error Pattern:**
```
<FACTOR_11> - Problème : la matrice est singulière ou presque singulière
```

**Fix Required:** Review test case or solver tolerance settings

---

## Priority Task Checklist

### High Priority (Dependency Issues)
- [x] ✅ **Add `medcoupling` to run dependencies** - Already fixed
- [ ] 🔴 **Fix MGIS/MFront configuration** - Set `ASTER_BEHAVIOUR_LIB` properly (20 tests)
- [ ] 🟡 Add `homard` to dependencies or mark tests as optional (17 tests)
- [ ] 🟡 Add `xmgrace` to dependencies or mark visualization tests as optional (5 tests)
- [ ] 🟡 Add `aspell` to dependencies (1 test)

### Medium Priority (Build/Configuration)
- [ ] 🟠 Fix Fortran compiler path for runtime compilation tests (1 test)
- [ ] 🟠 Configure PETSc or mark tests as expected failures (3 tests)
- [ ] 🟠 Fix MISS3D permissions or mark as expected failure (1 test)
- [ ] 🟠 Add or remove asrun dependency (2 tests)

### Low Priority (Debugging Required)
- [ ] 🔵 Debug segmentation faults (6 tests) - requires detailed analysis
- [ ] 🔵 Investigate numerical differences in NOOK_TEST_RESU tests (2 tests)
- [ ] 🔵 Review singular matrix issue (1 test)

---

## Recommendations

### Immediate Actions
1. **Fix MGIS Configuration**: Set `ASTER_BEHAVIOUR_LIB` in the configuration file or build process
2. **Add Optional Dependencies**: Create a conda metapackage with optional tools (homard, xmgrace, aspell)
3. **Document Known Limitations**: Clearly document which features are not available (PETSc, MISS3D)

### Short-term Actions
4. **Fix Compiler Paths**: Ensure runtime compilation uses available compilers, not build-time paths
5. **Test Categorization**: Mark tests requiring optional dependencies with proper pytest markers
6. **CI/CD Integration**: Set up automated tracking of test pass rates

### Long-term Actions
7. **Debug Segfaults**: Use valgrind/gdb to identify memory corruption issues
8. **Enable PETSc**: Consider adding PETSc solver support to the conda build
9. **Comprehensive Testing**: Add more granular test categories for better failure tracking

---

## Estimated Impact of Fixes

| Fix | Tests Affected | Expected Pass Rate After Fix |
|-----|----------------|------------------------------|
| Current | - | 96.18% |
| + MGIS config fix | 20 | 97.06% |
| + Add homard | 17 | 97.85% |
| + Add xmgrace | 5 | 98.07% |
| + Add aspell | 1 | 98.11% |
| + Fix segfaults | 6 | 98.37% |
| **Total Potential** | **49 tests** | **~98.4%** |

---

## Notes

- The medcoupling issue has already been addressed by adding it to run dependencies
- Some failures (PETSc, MISS3D, HOMARD) may be expected if those features are not intended for the conda distribution
- The 6 segmentation faults require individual investigation and may indicate deeper issues
- Overall, the 96.18% pass rate is very good for an initial conda packaging effort

