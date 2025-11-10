# Debugging Strategy for Empty Matrix Issue

## Quick Win: Add Debug Output

Modify `bibfor/third_party_interf/amumpm.F90` around line 510 to add extensive debug output:

```fortran
! After line 508 where nz2 = to_mumps_int(nzloc)
! Add this debug block:

write (ifm, *) '=== AMUMPM DEBUG START ==='
write (ifm, *) 'Matrix dimension n:', n
write (ifm, *) 'Total matrix terms (kterm):', kterm  
write (ifm, *) 'Non-zeros found (nzloc):', nzloc
write (ifm, *) 'Filtered (exactly zero):', nfilt1
write (ifm, *) 'Filtered (underflow):', nfilt3
write (ifm, *) 'Filtered (overflow):', nfilt2
write (ifm, *) 'Filter threshold (rfiltr):', rfiltr
write (ifm, *) 'rmin:', rmin
write (ifm, *) 'rmax:', rmax
write (ifm, *) 'epsmat:', epsmat
write (ifm, *) '=== AMUMPM DEBUG END ==='

! Then continue with the existing check:
maxnz2 = nzloc
```

This will tell you:
- Is the matrix completely empty? (kterm = 0)
- Are all values zero? (nfilt1 = kterm, nzloc = 0)
- Are values being filtered? (nfilt2 or nfilt3 large)
- Is the filter threshold wrong? (rfiltr too large)

## Compile and Test

```cmd
# Rebuild just the affected file
pixi run buildd --no-clean

# Install
pixi run installd

# Run a failing test
run_ctest --resutest=temp\debug -R ssnp167f --verbose

# Check output
type temp\debug\ssnp167f.mess | findstr "AMUMPM DEBUG"
```

## Interpret Results

### Case 1: kterm = 0
**Problem**: Matrix isn't being assembled at all
**Next**: Check matrix assembly routines before `amumpm` is called

### Case 2: nfilt1 = kterm (all zeros)
**Problem**: Matrix values are all zero
**Next**: Check element stiffness matrix computation

### Case 3: nfilt2 large (overflow)
**Problem**: Matrix values are corrupted (NaN, Inf, or huge numbers)
**Next**: Check for memory corruption or uninitialized values

### Case 4: nfilt3 large (underflow)
**Problem**: Matrix values are valid but too small
**Next**: Check scaling or check if rfiltr is too aggressive

### Case 5: nzloc > 0 but still fails
**Problem**: Issue is actually in MUMPS, not matrix assembly
**Next**: Check MUMPS integer size and calling convention

## Alternative: Test Simpler Case

Try a much simpler test that should definitely have a non-empty matrix:

```cmd
# Run a simple linear static test
run_ctest --resutest=temp\debug -R ssll100a --verbose
```

If even simple linear tests fail with FACTOR_41, it's a fundamental issue.
If they pass, it's specific to non-linear analysis.

## Check MUMPS Integer Type Mismatch

The `to_mumps_int` function converts Code_Aster integers to MUMPS integers.
If there's a size mismatch, this could corrupt data.

Check in your MUMPS build:
```bat
# In C:\AibelProgs\code\condapackaging\src\code_aster\mumps
# Look for Makefile.inc or build.bat

# Should have:
OPTF = /4I8 /4R8     # 8-byte integers and reals
OPTC = -DINTSIZE64   # or -DINTSIZE=64
```

If MUMPS was compiled with 4-byte integers (INT32) but Code_Aster expects 8-byte (INT64), this would cause:
- Index corruption
- Matrix data corruption
- Empty matrices being detected

## Check Calling Convention

Intel Fortran on Windows uses:
- `/names:lowercase` - Lowercase function names
- `/assume:underscore` - Add trailing underscore

Your MUMPS **must** be compiled with the same settings, or function calls will fail silently.

Test this:
```cmd
# Check what symbols MUMPS exports
dumpbin /EXPORTS .pixi\envs\debug\Library\bin\dmumps.dll | findstr "dmumps"

# Should see:
#   dmumps_    (lowercase with underscore)
# NOT:
#   DMUMPS_    (uppercase)
#   dmumps     (no underscore)
```

## Worst Case: Memory Corruption from JEVEUX

If removing `/check:bounds` masked a real memory corruption bug, the matrix might be getting corrupted during assembly.

Test by temporarily re-enabling `/check:bounds` for just one test:

```cmd
# In test command file, add at the top:
import os
os.environ['FCFLAGS'] = os.environ.get('FCFLAGS', '') + ' /check:bounds'

# Then run the test - if it fails with subscript error before reaching MUMPS,
# confirms memory corruption
```

## Summary of Debugging Path

1. ✅ **Add debug output** → Identify which case (1-5 above)
2. ✅ **Test simple case** → Confirm it's widespread or specific
3. ✅ **Check MUMPS build** → Verify INT64 and calling convention
4. ✅ **Check symbols** → Confirm MUMPS exports are correct
5. ⚠️ **Consider memory corruption** → Re-enable checks selectively

Once you know which case applies, the fix will be much clearer!

