@echo off
setlocal
set CLICOLOR_FORCE=1

call conda_env.bat

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

set LDFLAGS=%LDFLAGS% /LIBPATH:%LIB_PATH_ROOT%/lib pthread.lib

REM link with MSVCRT.lib
set CFLAGS=%CFLAGS% /MD

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
  --disable-petsc ^
  --maths-libs=auto ^
  --without-hg

waf install_debug -v

endlocal