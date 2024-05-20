@echo off
setlocal

set "MINIF_VERSION=24.3.0-0"
set "MINIF_EXE=Miniforge3-%MINIF_VERSION%-Windows-x86_64.exe"

curl -o %MINIF_EXE% -L https://github.com/conda-forge/miniforge/releases/download/%MINIF_VERSION%/%MINIF_EXE% && start /wait "" %MINIF_EXE% /InstallationType=JustMe /RegisterPython=0 /S /D=C:\Work\Miniforge3

endlocal