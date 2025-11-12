@echo off

set CLICOLOR_FORCE=1

setlocal enabledelayedexpansion

SET CONFIG_PARAMETERS_addmem=5000

echo "Setting compiler env vars"

:: set FC=flang-new.exe
set FC=ifx.exe
set CC=clang-cl.exe
set CXX=clang-cl.exe

SET OUTPUT_DIR=%SRC_DIR%/build
echo "OUTPUT_DIR: %OUTPUT_DIR%"

set ASTER_PLATFORM_MSVC64=1
set ASTER_PLATFORM_WINDOWS=1

set MKLROOT=%LIBRARY_PREFIX%
SET MKLROOT=%MKLROOT:\=/%

SET LIB_PATH_ROOT=%LIBRARY_PREFIX:\=/%
SET PREF_ROOT=%PREFIX:\=/%

set LIBPATH_HDF5=%LIB_PATH_ROOT%/lib
set INCLUDES_HDF5=%LIB_PATH_ROOT%/include

set LIBPATH_MED=%LIB_PATH_ROOT%/lib
set INCLUDES_MED=%LIB_PATH_ROOT%/include

set LIBPATH_METIS=%LIB_PATH_ROOT%/lib
set INCLUDES_METIS=%LIB_PATH_ROOT%/include

set LIBPATH_MUMPS=%LIB_PATH_ROOT%/lib
set "INCLUDES_MUMPS=%LIB_PATH_ROOT%/include %LIB_PATH_ROOT%/include/mumps_seq"

set LIBPATH_SCOTCH=%LIB_PATH_ROOT%/lib
set INCLUDES_SCOTCH=%LIB_PATH_ROOT%/include

set TFELHOME=%LIB_PATH_ROOT%

set LIBPATH_MGIS=%LIB_PATH_ROOT%/bin
set INCLUDES_MGIS=%LIB_PATH_ROOT%/include

REM Compiler flags
set LIBPATH=%PREF_ROOT%/libs %LIBPATH%

REM /MD link with MSVCRT.lib. /FS allow for c compiler calls to vc140.pdb on multiple threads (for cl.exe only)

set CFLAGS=%CFLAGS% /FS /MD
set CXXFLAGS=%CXXFLAGS% /MD

echo "Using Intel Fortran LLVM IFX compiler"
set FC_SEARCH=ifort
set FCFLAGS=%FCFLAGS% /fpp /integer-size:64 /real-size:64 /MD /names:lowercase /assume:underscore /assume:nobscc /fpe:0 /traceback /nologo
:: Add lib paths
set LDFLAGS=%LDFLAGS% /LIBPATH:%LIB_PATH_ROOT%/lib /LIBPATH:%LIB_PATH_ROOT%/bin /LIBPATH:%PREF_ROOT%/libs
:: Set up paths for Intel Fortran compiler in conda environment
set "PATH=%BUILD_PREFIX%\Library\bin;%BUILD_PREFIX%\Scripts;%PATH%"
set "LIB=%BUILD_PREFIX%\Library\lib;%LIB%"
set "INCLUDE=%BUILD_PREFIX%\opt\compiler\include\intel64;%BUILD_PREFIX%\Library\include;%INCLUDE%"

:: Signal to ifort.py that we're using conda-based Intel Fortran
set "INTEL_FORTRAN_VERSION=2025.1162"
set "CONDA_BUILD_INTEL_FORTRAN=1"

:: Increase compiler memory limits to avoid out-of-memory errors
@REM set "FOR_STACK_LIMIT=1000000000"
@REM set "_INTEL_COMPILER_HEAP_SIZE=2048"

if %CC% == "cl.exe" set CFLAGS=%CFLAGS% /sourceDependencies %OUTPUT_DIR%

:: Create dll debug pdb
if "%build_type%" == "debug" (
    echo "Building debug version"
    set FCFLAGS=%FCFLAGS% /check:stack /Z7 /traceback
    set CFLAGS=%CFLAGS% /Z7 /FS /Oy-
    set CXXFLAGS=%CXXFLAGS% /Z7 /FS /Oy-
    set LDFLAGS=%LDFLAGS% /DEBUG:FULL /INCREMENTAL:NO /OPT:REF /OPT:ICF /PDBALTPATH:%%_PDB%%
) else (
    echo "Building release version"
)

:: Add Math libs
set LDFLAGS=%LDFLAGS% mkl_intel_lp64_dll.lib mkl_intel_thread_dll.lib mkl_core_dll.lib libiomp5md.lib

:: Add threading libs
@REM set LDFLAGS=%LDFLAGS% pthread.lib

:: Add mumps libs
set LDFLAGS=%LDFLAGS% mpiseq.lib esmumps.lib scotch.lib scotcherr.lib scotcherrexit.lib

:: Add metis libs
set LDFLAGS=%LDFLAGS% metis.lib

:: Add libmed libs
set LDFLAGS=%LDFLAGS% med.lib medC.lib medfwrap.lib medimport.lib

set INCLUDES_BIBC=%PREF_ROOT%/include %SRC_DIR%/bibfor/include %INCLUDES_BIBC%

set DEFINES=H5_BUILT_AS_DYNAMIC_LIB _CRT_SECURE_NO_WARNINGS _SCL_SECURE_NO_WARNINGS WIN32_LEAN_AND_MEAN ASTER_PLATFORM_MSVC64

set DEFINES=%DEFINES% ASTER_INT8
if "%int_type%" == "64" (
    echo "Using 64-bit integer type"
) else (
    echo "Using 32-bit integer type"
)

python %RECIPE_DIR%\config\update_version.py

REM Install for standard sequential
waf configure ^
  --safe ^
  --check-fortran-compiler=ifort ^
  --use-config-dir=%SRC_DIR%/config/ ^
  --med-libs="med medC medfwrap medimport" ^
  --prefix=%LIB_PATH_ROOT% ^
  --out="%SRC_DIR%/build" ^
  --libdir="%LIBRARY_PREFIX%/lib" ^
  --bindir="%LIBRARY_PREFIX%/bin" ^
  --spdir=%SP_DIR% ^
  --disable-aster-subdir ^
  --enable-med ^
  --enable-hdf5 ^
  --enable-mumps ^
  --enable-metis ^
  --enable-scotch ^
  --enable-mfront ^
  --disable-mpi ^
  --disable-openmp ^
  --disable-petsc ^
  --install-tests ^
  --maths-libs=auto ^
  --msvc-entry ^
  --without-hg ^
  --without-repo

if errorlevel 1 (
    type %SRC_DIR%/build/%build_type%/config.log
    exit 1
)

if "%build_type%" == "debug" (
    waf install_debug
) else (
    waf install
)

if errorlevel 1 exit 1

REM Copy PDB files to the package for debugging support
echo Copying PDB files for debugging...
if exist "%LIBRARY_PREFIX%\lib\aster\*.pdb" (
    copy "%LIBRARY_PREFIX%\lib\aster\*.pdb" "%LIBRARY_BIN%\" >nul 2>&1
)
if exist "%LIBRARY_PREFIX%\bin\*.pdb" (
    echo PDB files already in bin directory
)

REM Note: We use /Z7 flag which embeds debug info in .obj files rather than separate PDBs
REM This makes the debug info self-contained and doesn't require source files to be packaged
REM Debuggers can still show source if the user has the source tree, but it's not required for symbols

REM Move code_aster and run_aster directories (including subdirectories)
@REM move "%LIBRARY_PREFIX%\lib\aster\code_aster" "%SP_DIR%\code_aster"
@REM move "%LIBRARY_PREFIX%\lib\aster\run_aster" "%SP_DIR%\run_aster"
@REM
@REM REM Move all .pyd files to %SP_DIR%
@REM for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.pyd") do move "%%f" "%SP_DIR%"
@REM
@REM REM Move all dll/lib/pdb files to %BIN_DIR%
@REM for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.dll") do move "%%f" "%LIBRARY_BIN%"
@REM for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.lib") do move "%%f" "%LIBRARY_BIN%"
@REM for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.pdb") do move "%%f" "%LIBRARY_BIN%"

echo Files moved successfully.

endlocal