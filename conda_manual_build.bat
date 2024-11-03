@echo off

setlocal enabledelayedexpansion

SET CONFIG_PARAMETERS_addmem=5000
SET PARENT_DIR=%~dp0
SET PARENT_DIR=%PARENT_DIR:\=/%

set INCLUDE_TESTS=0
set USE_LOG=0
set COLOR_ENABLED=1
:: BUILD_TYPE can be either debug or release
set BUILD_TYPE=debug
set CLEAN_BUILD=1
set PIXI_BUILD=0

:parse_args
if "%~1"=="" goto end_parse_args
if /i "%~1"=="--pixi-build" set PIXI_BUILD=1
if /i "%~1"=="--install-tests" set INCLUDE_TESTS=1
if /i "%~1"=="--use-log" set USE_LOG=1
if /i "%~1"=="--no-color" set COLOR_ENABLED=0
if /i "%~1"=="--no-clean" set CLEAN_BUILD=0
shift
goto parse_args

:end_parse_args

if %COLOR_ENABLED%==1 (
    echo Enabling color output
    set CLICOLOR_FORCE=1
) else (
    echo Disabling color output
    set CLICOLOR_FORCE=0
)

echo "Setting compiler env vars"
set "CC=clang-cl.exe"
set "CXX=clang-cl.exe"
set "FC=ifx.exe"
@REM set "LINK=link.exe"

:: TO set the number of cores, use the env variable JOBS
if %PIXI_BUILD% == 1 (
    echo "Using pixi build"
    echo CONDA_PREFIX=%CONDA_PREFIX%
    set "PREFIX=%CONDA_PREFIX%"
    set "LIBRARY_PREFIX=%CONDA_PREFIX%/Library"
    set "ASTER_ROOT=%CONDA_PREFIX%"
    set "RUNASTER_ROOT=%CONDA_PREFIX%/Library"
    REM call %PARENT_DIR%\ifx_env.bat
) else (
    call %PARENT_DIR%\conda_env.bat
)

if errorlevel 1 exit 1

if defined JOBS (
    echo "Using %JOBS% cores"
)

where python
where "%CC%"
where "%FC%"
REM where "%LINK%"

SET OUTPUT_DIR=%PARENT_DIR%/build/std
SET OUTPUT_DIR=%OUTPUT_DIR:\=/%

set ASTER_PLATFORM_MSVC64=1
REM set ASTER_BLAS_INT_SIZE=8
set ASTER_PLATFORM_WINDOWS=1
set ASTER_HAVE_MPI=0

set MKLROOT=%LIBRARY_PREFIX%
SET MKLROOT=%MKLROOT:\=/%

SET LIB_PATH_ROOT=%LIBRARY_PREFIX:\=/%
SET PREF_ROOT=%PREFIX:\=/%

REM set INCLUDES_PYBIND11=%LIB_PATH_ROOT%/include

set LIBPATH_HDF5=%LIB_PATH_ROOT%/lib
set INCLUDES_HDF5=%LIB_PATH_ROOT%/include

set LIBPATH_MED=%LIB_PATH_ROOT%/lib
set INCLUDES_MED=%LIB_PATH_ROOT%/include

set LIBPATH_METIS=%LIB_PATH_ROOT%/lib
set INCLUDES_METIS=%LIB_PATH_ROOT%/include

set LIBPATH_MUMPS=%LIB_PATH_ROOT%/lib
set "INCLUDES_MUMPS=%LIB_PATH_ROOT%/include"

set LIBPATH_SCOTCH=%LIB_PATH_ROOT%/lib
set INCLUDES_SCOTCH=%LIB_PATH_ROOT%/include

set TFELHOME=%LIB_PATH_ROOT%

set LIBPATH_MGIS=%LIB_PATH_ROOT%/bin
set INCLUDES_MGIS=%LIB_PATH_ROOT%/include

REM Compiler flags
set LIBPATH=%PREF_ROOT%/libs %LIBPATH%

REM /MD link with MSVCRT.lib. /FS allow for c compiler calls to vc140.pdb on multiple threads (for cl.exe only)

set CFLAGS=%CFLAGS% /FS /MD -Wno-visibility
set CXXFLAGS=%CXXFLAGS% /MD

if "%FC%" == "ifx.exe" (
    echo "Using Intel Fortran LLVM IFX compiler"
    set FC_SEARCH=ifort
    set FCFLAGS=%FCFLAGS% /fpp /MD /4I8 /4R8 /names:lowercase /assume:underscore /assume:nobscc /fpe:0
    :: Add lib paths
    set LDFLAGS=%LDFLAGS% /LIBPATH:%LIB_PATH_ROOT%/lib /LIBPATH:%LIB_PATH_ROOT%/bin /LIBPATH:%PREF_ROOT%/libs

) else (
    echo "Using LLVM Flang Fortran compiler"
    set FC_SEARCH=flang
    set FCFLAGS=%FCFLAGS% -cpp --dependent-lib=msvcrt -fdefault-double-8 -fdefault-integer-8 -fdefault-real-8 -funderscoring
    :: Add lib paths
    set LDFLAGS=%LDFLAGS% -L %LIB_PATH_ROOT%/lib -L %LIB_PATH_ROOT%/bin -L %PREF_ROOT%/libs
)

:: Create dll debug pdb
if "%BUILD_TYPE%" == "debug" (
    set FCFLAGS=%FCFLAGS% /check:stack
    set CFLAGS=%CFLAGS% /Zi
    set CXXFLAGS=%CXXFLAGS% /Zi
) else (
    REM set the equivalent of RelWithDebInfo
@REM     set FCFLAGS=%FCFLAGS% /debug:full /debug-parameters:all /traceback
@REM     set CFLAGS=%CFLAGS% /Zi
@REM     set CXXFLAGS=%CXXFLAGS% /Zi
)

if %CC% == "cl.exe" set CFLAGS=%CFLAGS% /sourceDependencies %OUTPUT_DIR%

:: Add Math libs
set LDFLAGS=%LDFLAGS% mkl_intel_lp64_dll.lib mkl_intel_thread_dll.lib mkl_core_dll.lib libiomp5md.lib

:: Add threading libs
set LDFLAGS=%LDFLAGS% pthread.lib

:: Add mumps libs
set LDFLAGS=%LDFLAGS% mpiseq.lib esmumps.lib scotch.lib scotcherr.lib scotcherrexit.lib

:: Add metis libs
set LDFLAGS=%LDFLAGS% metis.lib

:: Add libmed libs
set LDFLAGS=%LDFLAGS% med.lib medC.lib medfwrap.lib medimport.lib

set INCLUDES_BIBC=%PREF_ROOT%/include %PARENT_DIR%/bibfor/include %INCLUDES_BIBC%

set DEFINES=H5_BUILT_AS_DYNAMIC_LIB _CRT_SECURE_NO_WARNINGS _SCL_SECURE_NO_WARNINGS WIN32_LEAN_AND_MEAN
if "%build_type%" == "debug" (
@REM     set DEFINES=%DEFINES% ASTER_DEBUG_ALL
)
REM Clean the build directory
if %CLEAN_BUILD%==1 (
    echo "Cleaning build directory"
    waf distclean
)

REM Update version
python conda\scripts\update_version.py

set BUILD=std
set EXTRA_ARGS=
if "%INCLUDE_TESTS%" == "1" (
    echo "Including tests"
    set "EXTRA_ARGS=--install-tests"
)

REM Install for standard sequential
if %CLEAN_BUILD%==1 (
    waf configure ^
      --python=%PYTHON% ^
      --check-fortran-compiler=%FC_SEARCH% ^
      --use-config-dir=%PARENT_DIR%/config/ ^
      --med-libs="med medC medfwrap medimport" ^
      --prefix=%LIB_PATH_ROOT% ^
      --out=%OUTPUT_DIR% ^
      --enable-med ^
      --enable-hdf5 ^
      --enable-mumps ^
      --enable-metis ^
      --enable-scotch ^
      --disable-mpi ^
      --disable-petsc ^
      --maths-libs=auto ^
      --msvc-entry ^
      --without-hg ^
      --without-repo %EXTRA_ARGS%
)
REM   --install-tests ^
if "%errorlevel%" == "1" (
    type "%OUTPUT_DIR%/config.log"
    exit 1
)

REM Conditional log handling
if %USE_LOG%==1 (
    set "datetimeString="
    call conda_datetime.bat
    waf install_debug -vvv > "install_debug_%datetimeString%.log" 2>&1
) else (
    if "%BUILD_TYPE%" == "debug" (
        waf install_debug -v
    ) else (
        waf install -v
    )
)

if errorlevel 1 exit 1

call conda_rearrange.bat

if errorlevel 1 exit 1

endlocal