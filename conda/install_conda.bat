@echo off

setlocal

set "MINIF_VERSION=24.3.0-0"
set "MINIF_EXE=Miniforge3-%MINIF_VERSION%-Windows-x86_64.exe"

SET PARENT_DIR=%~dp0
SET TMP_DIR=%PARENT_DIR%temp
SET LOCAL_EXE=%TMP_DIR%\%MINIF_EXE%

if not exist %LOCAL_EXE% (
    if not exist %TMP_DIR% mkdir %TMP_DIR%
    echo "Downloading Miniforge3"
    curl -o %LOCAL_EXE% -L https://github.com/conda-forge/miniforge/releases/download/%MINIF_VERSION%/%MINIF_EXE%
else (
    echo "Miniforge3 already downloaded"
)

echo "Installing Miniforge3"
start /wait "" %LOCAL_EXE% /InstallationType=JustMe /RegisterPython=0 /S /D=C:\Work\Miniforge3

endlocal