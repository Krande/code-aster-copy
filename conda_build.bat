@echo off
setlocal
set CLICOLOR_FORCE=1

rem Set the python library prefix
set PYTHON_ENV=codeaster-deps
set LIBRARY_PREFIX=C:\Work\mambaforge\envs\%PYTHON_ENV%\Library

REM Set the path to the VS Cl and Intel fortran compiler
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
  set "CC=cl.exe"
  set "CXX=cl.exe"
  set "FC=ifort.exe"
)
where python
where cl
where ifort

SET PARENT_DIR=%~dp0
echo PARENT_DIR=%PARENT_DIR%
rem convert to forward slashes
SET PARENT_DIR=%PARENT_DIR:\=/%

set ASTER_PLATFORM_MSVC=1

set MKLROOT=%LIBRARY_PREFIX%
SET MKLROOT=%MKLROOT:\=/%

waf distclean

REM Install for standard sequential
waf configure ^
  --use-config-dir=%PARENT_DIR%/conda/ ^
  --prefix=%LIBRARY_PREFIX% ^
  --libdir=%LIBRARY_PREFIX%\libs\python312.lib ^
  --pythondir=%LIBRARY_PREFIX% ^
  --disable-mpi ^
  --install-tests ^
  --maths-libs=auto ^
  --embed-metis ^
  --without-hg

waf install_debug

endlocal