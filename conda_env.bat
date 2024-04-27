rem "INTEL_VARS_PATH=C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env"
rem "VS_VARS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build"
rem "CONDA_ROOT=C:\Work\miniconda3"

REM set the local path variables, INTEL_VARS_PATH, VS_VARS_PATH and CONDA_ROOT from .env file
for /f "tokens=*" %%a in (.env) do set %%a

rem Set the python library prefix
set PREFIX=%CONDA_ROOT%\envs\%PYTHON_ENV_NAME%
set LIBRARY_PREFIX=%PREFIX%\Library

REM Activate python env, env variables for VS Cl (or clang-cl) and Intel fortran compiler
@call "%CONDA_ROOT%\Scripts\activate.bat" %PYTHON_ENV_NAME%
call "%VS_VARS_PATH%\vcvars64.bat"
@call "%INTEL_VARS_PATH%\vars.bat" -arch intel64 vs2022

REM if variable "print" is passed, call printenv
if "%1" == "print" (
    set
)