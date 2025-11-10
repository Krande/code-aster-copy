# Next Steps for Code_Aster MSVC Windows Port

## Current Status (2025-11-10)
- ✅ Fixed `waitstatus_to_exitcode` OverflowError in `run_aster/utils.py`
- ✅ 406/2275 tests passing (17.85%)
- ❌ 1,869 tests failing with access violations (82.15%)

## Immediate Investigation Steps

### Step 1: Regenerate Symbol Export Files

The crashes are likely due to missing or incorrect symbol exports. Regenerate the `.def` files:

```cmd
pixi run def-bibfor
pixi run def-bibcpp
```

Then rebuild and reinstall:

```cmd
pixi run buildd --no-clean
pixi run installd
```

### Step 2: Identify Root Cause with Minimal Test

Pick the simplest failing test and debug it:

```cmd
cd temp\debug

REM Look for a simple static analysis test
findstr /c:"MECA_STATIQUE" *.mess | findstr /v "WIN_FATAL"
```

Good candidates for debugging:
- `ssll100b` - Simple linear static test
- `ssls100a` - Simple linear statics
- `erreu01a` - Error handling test (might show interface issues)

### Step 3: Run Test with Detailed Logging

Enable verbose logging for a single test:

```cmd
REM Run single test with Python debugging
cd C:\AibelProgs\code\code-aster-src

REM Set environment for detailed crash info
set PYTHONFAULTHANDLER=1

REM Run a simple failing test
run_ctest --resutest=temp\debug --verbose -R ssll100b
```

### Step 4: Check for Missing Symbols

Use `dumpbin` to check what symbols are exported:

```cmd
REM Check what's exported from the main library
dumpbin /EXPORTS .pixi\envs\debug\Library\bin\aster.dll > symbol_check.txt

REM Check for Fortran symbols
findstr /i "meca_statique\|calcul\|op0071" symbol_check.txt
```

### Step 5: Verify Fortran/C Interface

Check if the redirect modules are correctly set up:

```cmd
REM List redirect modules
dir msvc\c_entrypoints\*.cpp

REM Check if they're being built
dir .pixi\envs\debug\Library\bin\op*.dll
```

## Debugging Strategy

### Option A: Use WinDbg (Recommended)

1. Install WinDbg from Windows SDK
2. Run test under debugger:

```cmd
windbg -g -G python -m run_aster.run_aster_main --test --comm=<test>.comm
```

3. When crash occurs, examine:
   - `k` - Call stack
   - `r` - Registers
   - `!analyze -v` - Automatic analysis

### Option B: Use Visual Studio Debugger

1. Open Visual Studio 2022
2. Debug → Attach to Process
3. Attach to Python process running test
4. When crash occurs, examine Call Stack window

### Option C: Enable More Detailed Logging

Add debugging to the Python/C interface:

```python
# Edit: .pixi\envs\debug\lib\site-packages\code_aster\Supervis\ExecuteCommand.py
# Around line 282 (_call_oper), add:

import sys
print(f"DEBUG: About to call {self.get_name()}", file=sys.stderr)
sys.stderr.flush()

result = self._call_oper(*args, **kwargs)

print(f"DEBUG: Returned from {self.get_name()}", file=sys.stderr)
sys.stderr.flush()
```

## Known Issues to Check

### Issue 1: Calling Convention Mismatch

Intel Fortran and MSVC may have different default calling conventions. Check:

```fortran
! In Fortran code, ensure explicit interface declarations
! Example: bibfor/op/op0071.F90

subroutine op0071() bind(c, name='op0071_')
    ! Force C calling convention
end subroutine
```

### Issue 2: Name Mangling

Fortran name mangling on Windows may differ:
- GCC/gfortran: `op0071_` (lowercase with underscore)
- Intel Fortran: Could be `OP0071`, `op0071`, `op0071_`, or `_op0071`

Check current mangling:

```cmd
dumpbin /SYMBOLS .pixi\envs\debug\Library\lib\aster.lib | findstr "op0071"
```

### Issue 3: Integer Size Mismatches

Even with INT64 build, check for:

```fortran
! Bad - implicit integer (might be 32-bit)
integer :: n

! Good - explicit 64-bit integer
integer(kind=8) :: n
```

### Issue 4: Stack Alignment

MSVC requires 16-byte stack alignment for SSE instructions. Check compiler flags:

```cmd
REM Should see /arch:SSE2 or similar
type build\int64\debug\compile_commands.json | findstr "arch"
```

## Testing Approach

### Quick Verification Tests

Run these to verify basic functionality:

```cmd
REM 1. Test Python imports
python -c "import code_aster; print(code_aster.__version__)"

REM 2. Test library loading
python -c "from code_aster import libaster; print('OK')"

REM 3. Test basic operation
python -c "from code_aster.Commands import DEBUT, FIN; DEBUT(); FIN()"
```

### Isolate Problem Area

Create a minimal test case:

```python
# test_minimal.py
from code_aster.Commands import *

DEBUT(CODE="OUI")

# Try basic mesh operations
mesh = LIRE_MAILLAGE(FORMAT='ASTER', UNITE=20)

# Try basic model definition
model = AFFE_MODELE(
    MAILLAGE=mesh,
    AFFE=_F(TOUT='OUI', PHENOMENE='MECANIQUE', MODELISATION='3D')
)

FIN()
```

Run it:

```cmd
python test_minimal.py
```

## Build Configuration to Check

### Verify Compiler Flags

Check if proper flags are set:

```cmd
type build\int64\debug\c4che\_cache.py | findstr "DEFINES\|CXXFLAGS\|FCFLAGS"
```

Should have:
- `/DWIN32` or `/D_WIN32`
- `/D_USE_MATH_DEFINES`
- `/bigobj` for large object files
- `/MD` or `/MDd` for runtime library

### Check Link Flags

```cmd
type build\int64\debug\c4che\_cache.py | findstr "LINKFLAGS"
```

Should have:
- `/INCREMENTAL:NO` for debug builds
- `/DEBUG` for debug symbols
- Proper library paths

## Expected Outcomes

After investigating these areas, you should find one of:

1. **Missing symbols** → Regenerate .def files
2. **Wrong calling convention** → Update function declarations
3. **Name mangling mismatch** → Fix symbol names in .def files
4. **Stack corruption** → Add stack protection flags
5. **Memory alignment** → Fix struct definitions

## Success Criteria

### Phase 1 (This Week)
- [ ] Identify exact crash location in one test
- [ ] Understand why the crash occurs
- [ ] Fix one category of failures

### Phase 2 (Next Week)
- [ ] Apply fix broadly
- [ ] Reduce failures by 25% (to ~1,400 failing)
- [ ] Document the root cause

### Phase 3 (This Month)
- [ ] Fix major categories (MODI_MAILLAGE, MECA_STATIQUE)
- [ ] Achieve 50% pass rate
- [ ] Create test suite for Fortran/C interface

## Resources

### Useful Commands

```cmd
REM View test results
type temp\debug\<test>.mess

REM Check build log
type build\int64\debug\config.log

REM List all DLLs
dir /s /b .pixi\envs\debug\Library\bin\*.dll

REM Check dependencies
dumpbin /DEPENDENTS .pixi\envs\debug\Library\bin\aster.dll
```

### Key Files

- `run_aster/utils.py` - Test execution (FIXED)
- `code_aster/Supervis/ExecuteCommand.py` - Command execution
- `msvc/c_entrypoints/*.cpp` - Fortran→C redirect modules
- `msvc/scripts/generate_def_*.py` - Symbol export generation
- `bibfor/` - Fortran source code
- `bibcxx/` - C++ interface code

### Documentation

- Intel Fortran Compiler docs: https://www.intel.com/content/www/us/en/docs/fortran-compiler/
- MSVC ABI: https://learn.microsoft.com/en-us/cpp/build/reference/
- Code_Aster developer guide: Check `doc/` directory

## Contact Points

If stuck, consider:
1. Reviewing Code_Aster GitLab issues for similar problems
2. Checking Intel Fortran + MSVC compatibility guides
3. Consulting Windows HPC community forums
4. Reviewing MinGW/MSYS2 Code_Aster port (if exists) for comparison

---

**Remember**: The 17.85% pass rate shows the build infrastructure works. The crashes are systematic, suggesting one or two root causes affect most tests. Fixing those could dramatically improve the pass rate.

