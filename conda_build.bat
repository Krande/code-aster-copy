@echo off
setlocal
set CLICOLOR_FORCE=1

rem Set the python library prefix
set PYTHON_ENV=codeaster-deps
set PREFIX=C:\Work\mambaforge\envs\%PYTHON_ENV%
set LIBRARY_PREFIX=%PREFIX%\Library

REM Set the path to the VS Cl (or clang-cl) and Intel fortran compiler
set "INTEL_VARS_PATH=C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env"
set "VS_VARS_PATH=C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build"

rem check if regular VS2022 exists if not check if build tools 2022 exists
if not exist "%VS_VARS_PATH%" (
  set "VS_VARS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build"
)

@call "C:\Work\mambaforge\Scripts\activate.bat" %PYTHON_ENV%
call "%VS_VARS_PATH%\vcvars64.bat"
@call "%INTEL_VARS_PATH%\vars.bat" -arch intel64 vs2022

rem if not exist CC
if not exist "%CC%" (
  echo "Setting compiler env vars"
  set "CC=clang-cl.exe"
  set "CXX=clang-cl.exe"
  set "FC=ifx.exe"
)

set FCFLAGS=%FCFLAGS% -fpp

where python
where cl
where ifort

SET PARENT_DIR=%~dp0
echo PARENT_DIR=%PARENT_DIR%
rem convert to forward slashes
SET PARENT_DIR=%PARENT_DIR:\=/%

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
set INCLUDES_MUMPS=%LIB_PATH_ROOT%/include

set LIBPATH_SCOTCH=%LIB_PATH_ROOT%/lib
set INCLUDES_SCOTCH=%LIB_PATH_ROOT%/include

set TFELHOME=%LIB_PATH_ROOT%

set LIBPATH_MGIS=%LIB_PATH_ROOT%/bin
set INCLUDES_MGIS=%LIB_PATH_ROOT%/include

set LINKFLAGS=%LINKFLAGS% /LIBPATH:%LIB_PATH_ROOT%/lib pthread.lib

set INCLUDES_BIBC=%PREF_ROOT%/include

set DEFINES=H5_BUILT_AS_DYNAMIC_LIB

waf distclean

REM Install for standard sequential
waf configure ^
  --use-config-dir=%PARENT_DIR%/conda/ ^
  --med-libs=medC ^
  --prefix=%LIBRARY_PREFIX% ^
  --disable-mpi ^
  --install-tests ^
  --maths-libs=auto ^
  --without-hg

waf install_debug -v

endlocal