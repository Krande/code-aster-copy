@echo off
setlocal

@call conda_env.bat

link.exe %*

endlocal