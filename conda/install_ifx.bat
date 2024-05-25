@echo off
setlocal

:: Install Intel Fortran Compiler 2024
SET PARENT_DIR=%~dp0
SET PARENT_DIR=%PARENT_DIR:\=/%
:: 2024.0
@REM SET INSTALL_EXE=w_fortran-compiler_p_2024.0.2.27_offline.exe
@REM SET INSTALL_VER=3a64aab4-3c35-40ba-bc9c-f80f136a8005
:: 2024.1
SET INSTALL_EXE=w_fortran-compiler_p_2024.1.0.466_offline.exe
SET INSTALL_VER=f6a44238-5cb6-4787-be83-2ef48bc70cba

if not exist %PARENT_DIR%\%INSTALL_EXE% (
    curl https://registrationcenter-download.intel.com/akdlm/IRC_NAS/%INSTALL_VER%/%INSTALL_EXE% -o %PARENT_DIR%\%INSTALL_EXE%
)

%PARENT_DIR%\%INSTALL_EXE% -s --a --eula=accept --log-dir=%PARENT_DIR%

endlocal