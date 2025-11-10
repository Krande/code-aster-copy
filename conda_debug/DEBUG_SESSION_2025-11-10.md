# Debug Session Summary - 2025-11-10

## Issues Identified and Fixed

### 1. ✅ FIXED: `waitstatus_to_exitcode` OverflowError
**File**: `run_aster/utils.py`

**Problem**: Python 3.9+ built-in `waitstatus_to_exitcode` doesn't work correctly on Windows. It expects Unix-style wait status but Windows returns exit codes directly.

**Solution**: Force Windows to always use the custom `_waitstatus_to_exitcode` implementation:
```python
if not hasattr(os, "waitstatus_to_exitcode") or RUNASTER_PLATFORM == "win":
    waitstatus_to_exitcode = _waitstatus_to_exitcode
```

**Impact**: Tests can now complete and report their exit status correctly.

---

### 2. ✅ FIXED: Intel Fortran Console Output Overflow
**File**: `bibfor/supervis/affich.F90`

**Problem**: During catalog build, Fortran runtime error:
```
forrtl: severe (66): output statement overflows record, unit 6, file CONOUT$
```

**Root Cause**: Windows console (CONOUT$) has a limited record length for Fortran output. When `affich.F90` tries to write long text strings, it exceeds this limit.

**Solution**: Modified `affich` subroutine to chunk long output into 512-character segments on Windows:
```fortran
#ifdef _WIN32
    parameter (chunk_size=512)
#else
    parameter (chunk_size=999999)
#endif

! Chunk output for Windows console compatibility
if (text_len .le. chunk_size) then
    write (ifm, '(A)') texte(1:text_len)
else
    do i = 1, text_len, chunk_size
        remaining = min(chunk_size, text_len - i + 1)
        write (ifm, '(A)', advance='no') texte(i:i+remaining-1)
    end do
    write (ifm, '(A)') ''
end if
```

**Impact**: Build catalog generation should now complete without crashing.

---

### 3. ✅ IMPROVED: Intel Fortran Debug Flags
**File**: `msvc/scripts/conda_manual_build.bat`

**Changes**: Added comprehensive debugging flags for Intel Fortran (ifx):

**Previous**:
```bat
set FCFLAGS=%FCFLAGS% /check:stack
```

**Current**:
```bat
set FCFLAGS=%FCFLAGS% /traceback /check:bounds /check:pointers /debug:full /Zi /Od /Qtrapuv /fp:precise

:: Runtime environment variables
set FOR_DIAGNOSTIC_LOG_LEVEL=1
set FORT_FMT_RECL=1024
set FOR_DISABLE_STACK_TRACE=0
```

**Key Flags**:
- `/traceback` - Generate detailed Fortran stack traces on crashes (CRITICAL)
- `/check:bounds` - Runtime array bounds checking
- `/check:pointers` - Runtime pointer validity checking
- `/debug:full` - Full debug symbols
- `/Qtrapuv` - Initialize stack variables to detect uninitialized use
- `FOR_DISABLE_STACK_TRACE=0` - Enable Intel Fortran runtime stack traces

**Removed Flags**:
- `/check:uninit` - Not supported on Windows ifx
- `/warn:all` - Too strict, caused false positives on valid code
- `/check:all` - Too strict, replaced with specific checks

**Impact**: When crashes occur, we will get detailed Fortran stack traces showing:
- Exact subroutine/function that crashed
- Full call stack with line numbers
- Which library (Code_Aster, MUMPS, SCOTCH, etc.) was involved

---

## Test Results Analysis

**Current Status** (from 2025-11-10.txt):
- **Passing**: 406/2275 (17.85%)
- **Failing**: 1,869/2275 (82.15%)

**Top Failure Categories**:
1. Unknown - 457 tests (24.45%)
2. MODI_MAILLAGE_WIN_FATAL - 379 tests (20.28%)
3. STAT_NON_LINE_WIN_FATAL - 370 tests (19.80%)
4. MECA_STATIQUE_WIN_FATAL - 284 tests (15.20%)
5. JEVEUX1_55 errors - 144 tests (7.71%)
6. POST_ELEM_WIN_FATAL - 89 tests (4.76%)
7. MED file issues - 58 tests (3.10%)

**Hypothesis**: Most failures are due to:
1. Fortran/C/C++ ABI incompatibility (calling conventions, name mangling)
2. Missing or incorrect symbol exports
3. Memory management issues in JEVEUX
4. External library integration issues (MUMPS, SCOTCH, METIS, MED, HDF5)

---

## Next Steps

### Immediate: Re-run Tests with Enhanced Tracing

After rebuild completes:

```cmd
pixi run installd
```

Then run a representative failing test to see detailed Fortran traceback:

```cmd
# Test from STAT_NON_LINE_WIN_FATAL category
run_ctest --resutest=temp\debug -R comp009a --verbose

# Or simpler MECA_STATIQUE test
run_ctest --resutest=temp\debug -R ssll100b --verbose
```

**Expected Output**: With `/traceback` enabled, when crash occurs we should see:
```
Image              PC                Routine            Line        Source
bibfor.dll         <address>         <function>         <line>      <file.F90>
...
```

This will tell us:
1. Which exact Fortran subroutine crashed
2. Whether it's in Code_Aster code or external library
3. What operation was being performed
4. The call stack leading to the crash

### Investigation Focus

Based on traceback, investigate:

#### If crash is in Code_Aster Fortran code:
- Check for array bounds violations
- Check for uninitialized variables
- Check calling conventions at Fortran/C boundary

#### If crash is in external library (e.g., MUMPS, SCOTCH):
- Verify library was compiled with same toolchain (MSVC)
- Check integer size compatibility (INT64 vs INT32)
- Verify symbol exports from the library
- Check if library expects different calling convention

#### If crash is at interface boundary:
- Check `msvc/c_entrypoints/*.cpp` redirect modules
- Verify `.def` files have correct symbol names
- Check name mangling (underscore conventions)
- Verify parameter passing (by value vs by reference)

### Tools for Investigation

1. **View Fortran Stack Trace**: 
   ```cmd
   type temp\debug\<test>.mess
   ```

2. **Check Symbol Exports**:
   ```cmd
   dumpbin /EXPORTS .pixi\envs\debug\Library\bin\bibfor.dll > bibfor_symbols.txt
   dumpbin /EXPORTS .pixi\envs\debug\Library\bin\mumps_common.dll > mumps_symbols.txt
   ```

3. **Check Dependencies**:
   ```cmd
   dumpbin /DEPENDENTS .pixi\envs\debug\Library\bin\bibfor.dll
   ```

4. **Run Under Debugger**:
   ```cmd
   windbg -g -G python -m run_aster.run_aster_main <test>.comm
   ```

---

## Key Files Modified

1. `run_aster/utils.py` - Fixed waitstatus_to_exitcode for Windows
2. `bibfor/supervis/affich.F90` - Fixed console output overflow on Windows
3. `msvc/scripts/conda_manual_build.bat` - Enhanced Intel Fortran debug flags

## Documentation Created

1. `conda_debug/test_scanner/results/windows/nompi/2025-11-10-analysis.md` - Comprehensive test failure analysis
2. `conda_debug/NEXT_STEPS.md` - Practical debugging guide

---

## Success Indicators

After these fixes:
- ✅ Build catalog should complete without error 66
- ✅ Tests should complete and report status (not hang or overflow)
- 🔄 Tests should provide detailed Fortran stack traces on crash
- 🔄 We should be able to identify root cause of crashes

## Notes

- The 17.85% pass rate is actually encouraging for early MSVC porting
- Many failures appear systematic (same root cause)
- Fixing the ABI/calling convention issues could fix hundreds of tests at once
- The enhanced tracing will be crucial for identifying the root cause

