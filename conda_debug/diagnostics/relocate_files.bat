@echo off

setlocal

REM Define source and destination directories
set "LIBRARY_PREFIX=%CONDA_PREFIX%\Library"
set "SP_DIR=%CONDA_PREFIX%\Lib\site-packages"
set "BIN_DIR=%LIBRARY_PREFIX%\bin"

REM Move code_aster and run_aster directories (including subdirectories)
move "%LIBRARY_PREFIX%\lib\aster\code_aster" "%SP_DIR%\code_aster"
move "%LIBRARY_PREFIX%\lib\aster\run_aster" "%SP_DIR%\run_aster"

REM Move all .pyd files to %SP_DIR%
for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.pyd") do move "%%f" "%SP_DIR%"

REM Move all dll/lib/pdb files to %BIN_DIR%
for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.dll") do move "%%f" "%BIN_DIR%"
for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.lib") do move "%%f" "%BIN_DIR%"
for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.pdb") do move "%%f" "%BIN_DIR%"

echo Files moved successfully.

endlocal