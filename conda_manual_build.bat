@echo off

setlocal enabledelayedexpansion

SET CONFIG_PARAMETERS_addmem=5000
SET PARENT_DIR=%~dp0
SET PARENT_DIR=%PARENT_DIR:\=/%

set INCLUDE_TESTS=0
set USE_LOG=0
set COLOR_ENABLED=1

:parse_args
if "%~1"=="" goto end_parse_args
if /i "%~1"=="--install-tests" set INCLUDE_TESTS=1
if /i "%~1"=="--use-log" set USE_LOG=1
if /i "%~1"=="--no-color" set COLOR_ENABLED=0
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
call %PARENT_DIR%\conda_env.bat

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
set ASTER_PLATFORM_WINDOWS=1
set ASTER_HAVE_MPI=0

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
set "INCLUDES_MUMPS=%LIB_PATH_ROOT%/include"

set LIBPATH_SCOTCH=%LIB_PATH_ROOT%/lib
set INCLUDES_SCOTCH=%LIB_PATH_ROOT%/include

set TFELHOME=%LIB_PATH_ROOT%

set LIBPATH_MGIS=%LIB_PATH_ROOT%/bin
set INCLUDES_MGIS=%LIB_PATH_ROOT%/include

REM Compiler flags
set LIBPATH=%PREF_ROOT%/libs %LIBPATH%

REM /MD link with MSVCRT.lib. /FS allow for c compiler calls to vc140.pdb on multiple threads (for cl.exe only)

set CFLAGS=%CFLAGS% /FS /MD /DMKL_ILP64
set CXXFLAGS=%CXXFLAGS% /MD /DMKL_ILP64

if "%FC%" == "ifx.exe" (
    echo "Using Intel Fortran LLVM IFX compiler"
    set FC_SEARCH=ifort
    set FCFLAGS=%FCFLAGS% /fpp /MD /4I8 /double-size:64 /real-size:64 /integer-size:64 /names:lowercase /assume:underscore /assume:nobscc
    :: Add lib paths
    set LDFLAGS=%LDFLAGS% /LIBPATH:%LIB_PATH_ROOT%/lib /LIBPATH:%LIB_PATH_ROOT%/bin /LIBPATH:%PREF_ROOT%/libs
) else (
    echo "Using LLVM Flang Fortran compiler"
    set FC_SEARCH=flang
    set FCFLAGS=%FCFLAGS% -cpp --dependent-lib=msvcrt -fdefault-double-8 -fdefault-integer-8 -fdefault-real-8 -funderscoring
    :: Add lib paths
    set LDFLAGS=%LDFLAGS% -L %LIB_PATH_ROOT%/lib -L %LIB_PATH_ROOT%/bin -L %PREF_ROOT%/libs
)
if %CC% == "cl.exe" set CFLAGS=%CFLAGS% /sourceDependencies %OUTPUT_DIR%

:: Create dll debug pdb
set LDFLAGS=%LDFLAGS% /DEBUG:FULL /INCREMENTAL:NO

:: Add Math libs
set LDFLAGS=%LDFLAGS% mkl_intel_ilp64_dll.lib mkl_intel_thread_dll.lib mkl_core_dll.lib libiomp5md.lib

:: Add threading libs
set LDFLAGS=%LDFLAGS% pthread.lib

:: Add mumps libs
set LDFLAGS=%LDFLAGS% mpiseq.lib esmumps.lib scotch.lib scotcherr.lib scotcherrexit.lib

:: Add metis libs
set LDFLAGS=%LDFLAGS% metis.lib

:: Add libmed libs
set LDFLAGS=%LDFLAGS% med.lib medC.lib medfwrap.lib medimport.lib

set INCLUDES_BIBC=%PREF_ROOT%/include %PARENT_DIR%/bibfor/include %INCLUDES_BIBC%

set DEFINES=H5_BUILT_AS_DYNAMIC_LIB PYBIND11_NO_ASSERT_GIL_HELD_INCREF_DECREF
REM Clean the build directory
@REM waf distclean

python conda\scripts\update_version.py

set BUILD=std

REM Install for standard sequential
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
  --enable-openmp ^
  --enable-metis ^
  --enable-scotch ^
  --disable-mpi ^
  --disable-petsc ^
  --install-tests ^
  --maths-libs=auto ^
  --without-hg ^
  --without-repo

REM Conditional log handling
if %USE_LOG%==1 (
    set "datetimeString="
    call conda_datetime.bat
    waf install_debug -vvv > "install_debug_%datetimeString%.log" 2>&1
) else (
    waf install_debug -v
)

endlocal