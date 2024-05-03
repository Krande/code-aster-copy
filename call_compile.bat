@echo off
setlocal

set PARENT_DIR=%~dp0

call %PARENT_DIR%\conda_env.bat

%*

endlocal