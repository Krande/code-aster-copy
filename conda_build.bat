@echo off
setlocal

set USE_LOG=0
set CLICOLOR_FORCE=1

call conda_env.bat

echo "Setting compiler env vars"
set "CC=cl.exe"
set "CXX=cl.exe"
set "FC=ifx.exe"
REM set "LINK_CC=XILINK.exe"
REM set "LINK_CXX=XILINK.exe"
REM set "LINK_FC=XILINK.exe"
REM set "AR=XILIB.exe"
REM set "RANLIB=XILIB.exe"
REM set "LD=XILINK.exe"

where python
where cl
where ifx

SET PARENT_DIR=%~dp0
SET PARENT_DIR=%PARENT_DIR:\=/%

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
set "INCLUDES_MUMPS=%LIB_PATH_ROOT%/include %LIB_PATH_ROOT%/include/mumps_seq"

set LIBPATH_SCOTCH=%LIB_PATH_ROOT%/lib
set INCLUDES_SCOTCH=%LIB_PATH_ROOT%/include

set TFELHOME=%LIB_PATH_ROOT%

set LIBPATH_MGIS=%LIB_PATH_ROOT%/bin
set INCLUDES_MGIS=%LIB_PATH_ROOT%/include

REM Compiler flags

REM /MD link with MSVCRT.lib. /FS allow for c compiler calls to vc140.pdb on multiple threads (for cl.exe only)
set CFLAGS=%CFLAGS% /FS /LD /sourceDependencies %OUTPUT_DIR%

set FCFLAGS=%FCFLAGS% -fpp /MD
set FCFLAGS=%FCFLAGS% /names:lowercase

set LDFLAGS=%LDFLAGS% /LIBPATH:%LIB_PATH_ROOT%/lib pthread.lib /DEBUG

set INCLUDES_BIBC=%PREF_ROOT%/include %PARENT_DIR%/bibfor/include %INCLUDES_BIBC%

set DEFINES=H5_BUILT_AS_DYNAMIC_LIB
REM Clean the build directory
waf distclean

REM set FORCE_BIBFOR_SEQUENCE=1
REM set MANUALLY_ADD_BIBFOR_DEPS=1

python conda\update_version.py

set BUILD=std

REM Install for standard sequential
waf configure ^
  --safe ^
  --check-fortran-compiler=ifort ^
  --use-config-dir=%PARENT_DIR%/config/ ^
  --med-libs=medC ^
  --prefix=%LIBRARY_PREFIX% ^
  --out=%OUTPUT_DIR% ^
  --disable-mpi ^
  --install-tests ^
  --maths-libs=auto ^
  --without-hg

REM if USE_LOG is set, then log the output to a file
if %USE_LOG%==1 (
    REM set a datetime variable down to the minute
    @call conda_datetime.bat
    waf install_debug -v > install_debug_%datetimeString%.log 2>&1
) else (
    waf install_debug -v -j 24
)

endlocal