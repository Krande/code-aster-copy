@echo off

setlocal enabledelayedexpansion

SET CONFIG_PARAMETERS_addmem=5000
SET PARENT_DIR=%~dp0../..
SET PARENT_DIR=%PARENT_DIR:\=/%

set INCLUDE_TESTS=0
set USE_LOG=0
set COLOR_ENABLED=1
:: BUILD_TYPE can be either debug or release
set BUILD_DEBUG=0
set BUILD64=0
set CLEAN_BUILD=1
set PIXI_BUILD=0
set VERBOSE_WAF=0
set EXTRA_WAF_ARGS=

:parse_args
if "%~1"=="" goto end_parse_args
if /i "%~1"=="--pixi-build" (
    set PIXI_BUILD=1
) else if /i "%~1"=="--install-tests" (
    set INCLUDE_TESTS=1
) else if /i "%~1"=="--use-log" (
    set USE_LOG=1
) else if /i "%~1"=="--no-color" (
    set COLOR_ENABLED=0
) else if /i "%~1"=="--no-clean" (
    set CLEAN_BUILD=0
) else if /i "%~1"=="--verbose-waf" (
    set VERBOSE_WAF=1
) else if /i "%~1"=="--int-type64" (
    set BUILD64=1
) else if /i "%~1"=="--build-debug" (
    set BUILD_DEBUG=1
) else (
    REM Collect unrecognized arguments
    set EXTRA_WAF_ARGS=!EXTRA_WAF_ARGS! %~1
)
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
    call %~dp0\conda_env.bat
)

if errorlevel 1 exit 1

if defined JOBS (
    echo "Using %JOBS% cores"
)

where python
where "%CC%"
where "%FC%"
REM where "%LINK%"

if %BUILD64% == 1 (
    set INTPATH=int64
) else (
    set INTPATH=int32
)

SET OUTPUT_DIR=%PARENT_DIR%/build/%INTPATH%

SET OUTPUT_DIR=%OUTPUT_DIR:\=/%

set ASTER_PLATFORM_MSVC64=1
REM set ASTER_BLAS_INT_SIZE=8
set ASTER_PLATFORM_WINDOWS=1
set ASTER_HAVE_MPI=0

set MKLROOT=%LIBRARY_PREFIX%
SET MKLROOT=%MKLROOT:\=/%

SET LIB_PATH_ROOT=%LIBRARY_PREFIX:\=/%
SET PREF_ROOT=%PREFIX:\=/%

set INCLUDES_PYBIND11=%LIB_PATH_ROOT%/include

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
set CXXFLAGS=%CXXFLAGS% /MD -Wno-visibility

if "%FC%" == "ifx.exe" (
    echo "Using Intel Fortran LLVM IFX compiler"
    set FC_SEARCH=ifort
    set FCFLAGS=%FCFLAGS% /4R8 /fpp /MD /names:lowercase /assume:underscore /assume:nobscc /fpe:0
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
if "%BUILD_DEBUG%" == "1" (
    echo "Setting debug flags"
    set FCFLAGS=%FCFLAGS% /check:stack
    set CFLAGS=%CFLAGS% /Zi
    set CXXFLAGS=%CXXFLAGS% /Zi
) else (
    echo "Setting release flags"
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

if %BUILD64% == 1 (
    echo "Using 64-bit integer type"
    if "%FC%" == "ifx.exe" (
        set FCFLAGS=%FCFLAGS% /4I8
    ) else (
        set FCFLAGS=%FCFLAGS% -fdefault-integer-8
    )
    set DEFINES=%DEFINES% ASTER_INT8
) else (
    echo "Using 32-bit integer type"
    if "%FC%" == "ifx.exe" (
        set FCFLAGS=%FCFLAGS% /4I8
    ) else (
        set FCFLAGS=%FCFLAGS% -fdefault-integer-8
    )
    set DEFINES=%DEFINES% ASTER_INT8
)

REM Clean the build directory
if %CLEAN_BUILD%==1 (
    echo "Cleaning build directory"
    waf distclean
)

REM Update version
python conda\scripts\update_version.py

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
      --enable-mfront ^
      --enable-scotch ^
      --disable-openmp ^
      --disable-mpi ^
      --disable-petsc ^
      --maths-libs=auto ^
      --msvc-entry ^
      --without-hg ^
      --without-repo %EXTRA_ARGS% %EXTRA_WAF_ARGS%
)
REM   --install-tests ^

if errorlevel 1 (
    type %OUTPUT_DIR%/config.log
    exit 1
)
set WAF_INSTALL_EXTRA_ARGS=
if %VERBOSE_WAF%==1 (
    set WAF_INSTALL_EXTRA_ARGS=-vvv
)
REM Conditional log handling
if %USE_LOG%==1 (
    set "datetimeString="
    call conda_datetime.bat
    waf install_debug -vvv > "install_debug_%datetimeString%.log" 2>&1
) else (
    if "%BUILD_DEBUG%" == "1" (
        waf install_debug %WAF_INSTALL_EXTRA_ARGS%
    ) else (
        waf install %WAF_INSTALL_EXTRA_ARGS%
    )
)

if errorlevel 1 exit 1

@REM call %~dp0\conda_rearrange.bat

if errorlevel 1 exit 1

endlocal