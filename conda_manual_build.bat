@echo off

setlocal enabledelayedexpansion

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


:: TO set the number of cores, use the env variable JOBS
call %PARENT_DIR%\conda_env.bat

if defined JOBS (
    echo "Using %JOBS% cores"
)

echo "Setting compiler env vars"
set "CC=clang-cl.exe"
set "CXX=clang-cl.exe"
set "FC=ifx.exe"
@REM set "LINK=link.exe"

where python
where "%CC%"
where "%FC%"
REM where "%LINK%"


SET OUTPUT_DIR=%PARENT_DIR%/build/std
SET OUTPUT_DIR=%OUTPUT_DIR:\=/%

set ASTER_PLATFORM_MSVC=1
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
set "INCLUDES_MUMPS=%LIB_PATH_ROOT%/include"

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
set FCFLAGS=%FCFLAGS% /fpp /MD
set FCFLAGS=%FCFLAGS% /names:lowercase /assume:underscore /assume:nobscc

if %CC% == "cl.exe" set CFLAGS=%CFLAGS% /sourceDependencies %OUTPUT_DIR%

:: Add lib paths
set LDFLAGS=%LDFLAGS% /LIBPATH:%LIB_PATH_ROOT%/lib /LIBPATH:%LIB_PATH_ROOT%/bin /LIBPATH:%PREF_ROOT%/libs

:: Add Math libs
set LDFLAGS=%LDFLAGS% mkl_intel_lp64_dll.lib mkl_intel_thread_dll.lib mkl_core_dll.lib libiomp5md.lib

:: Add windows libs
set LDFLAGS=%LDFLAGS% kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib

:: Add threading libs
set LDFLAGS=%LDFLAGS% pthread.lib

:: Add hdf5 libs
@REM set LDFLAGS=%LDFLAGS% hdf5.lib hdf5_hl.lib

:: Add mumps libs
set LDFLAGS=%LDFLAGS% mpiseq.lib scotch.lib scotcherr.lib metis.lib pord.lib esmumps.lib

:: Add libmed libs
set LDFLAGS=%LDFLAGS% med.lib medC.lib medfwrap.lib medimport.lib

set INCLUDES_BIBC=%PREF_ROOT%/include %PARENT_DIR%/bibfor/include %INCLUDES_BIBC%

set DEFINES=H5_BUILT_AS_DYNAMIC_LIB PYBIND11_NO_ASSERT_GIL_HELD_INCREF_DECREF
REM Clean the build directory
waf distclean

python conda\update_version.py

set BUILD=std

REM Install for standard sequential
waf configure ^
  --safe ^
  --check-fortran-compiler=ifort ^
  --use-config-dir=%PARENT_DIR%/config/ ^
  --med-libs="med medC medfwrap medimport" ^
  --prefix=%LIB_PATH_ROOT% ^
  --out=%OUTPUT_DIR% ^
  --embed-mumps ^
  --disable-mpi ^
  --maths-libs=auto ^
  --install-tests ^
  --without-hg


@REM   --disable-openmp ^
@REM   --disable-mumps ^

REM Conditional log handling
if %USE_LOG%==1 (
    set "datetimeString="
    call conda_datetime.bat
    waf install_debug -vvv > "install_debug_%datetimeString%.log" 2>&1
) else (
    waf install_debug -v
)

endlocal