@echo off
setlocal

set PARENT_DIR=%~dp0

call %PARENT_DIR%\ifx_env.bat

%*

endlocal