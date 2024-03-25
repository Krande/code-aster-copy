@echo off
setlocal

set USE_LOG=0
set CLICOLOR_FORCE=1

call conda_env.bat

link.exe %*

endlocal