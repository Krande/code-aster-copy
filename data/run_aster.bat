echo OFF
setlocal
set PYTHONIOENCODING=UTF-8
chcp 65001

:: If not CONDA_PREFIX is defined, set these variables
if not defined CONDA_PREFIX (
    set "RUNASTER_ROOT=%~dp0.."
    set "ASTER_ROOT=%RUNASTER_ROOT%\.."
    set "OUTILS=%ASTER_ROOT%\outils"
    set "PYTHONHOME=%ASTER_ROOT%"
    set "PYTHONPATH=%ASTER_ROOT%\lib\site-packages;%RUNASTER_ROOT%\lib\aster"
    set "PATH=%PYTHONHOME%;%OUTILS%;%PATH%"
) else (
    set "ASTER_ROOT=%CONDA_PREFIX%\Library"

)

call "%ASTER_ROOT%\share\aster\profile.bat

python -m run_aster.run_aster_main %*
