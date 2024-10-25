@echo off

rem "INTEL_VARS_PATH=C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env"
rem "VS_VARS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build"
rem "CONDA_ROOT=C:\Work\miniconda3"
SET PARENT_DIR=%~dp0
SET PARENT_DIR=%PARENT_DIR:\=/%
REM set the local path variables, INTEL_VARS_PATH, VS_VARS_PATH and CONDA_ROOT from .env file
for /f "tokens=*" %%a in (%PARENT_DIR%/.env) do set %%a

REM Activate python env, env variables for VS Cl (or clang-cl) and Intel fortran compiler
REM if env var DONOT_ACTIVATE is set, we can just exit from this batch script now

call "%VS_VARS_PATH%\vcvars64.bat"
if "%FC%" == "ifx.exe" (
    echo "Activating Intel Fortran compiler"
    echo "Vars bat: %INTEL_VARS_PATH%\vars.bat"
    @call "%INTEL_VARS_PATH%\vars.bat" -arch intel64
    set "CC=clang-cl.exe"
    set "CXX=clang-cl.exe"
    set "FC=%INTEL_VARS_PATH%/bin/ifx.exe"
)

REM if variable "print" is passed, call printenv
if "%1" == "print" (
    set
)