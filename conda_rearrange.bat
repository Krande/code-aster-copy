setlocal enabledelayedexpansion

if not defined SP_DIR (
    set SP_DIR=%CONDA_PREFIX%\Lib\site-packages
)
if not defined LIBRARY_PREFIX (
    set LIBRARY_PREFIX=%CONDA_PREFIX%\Library
)
if not defined LIBRARY_BIN (
    set LIBRARY_BIN=%LIBRARY_PREFIX%\bin
)

REM Move code_aster and run_aster directories (including subdirectories)
move "%LIBRARY_PREFIX%\lib\aster\code_aster" "%SP_DIR%\code_aster"
move "%LIBRARY_PREFIX%\lib\aster\run_aster" "%SP_DIR%\run_aster"

REM Move all .pyd files to %SP_DIR%
for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.pyd") do move "%%f" "%SP_DIR%"

REM Move all dll/lib/pdb files to %BIN_DIR%
for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.dll") do move "%%f" "%LIBRARY_BIN%"
for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.lib") do move "%%f" "%LIBRARY_BIN%"
for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.pdb") do move "%%f" "%LIBRARY_BIN%"

echo Files moved successfully.

endlocal