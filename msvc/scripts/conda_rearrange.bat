@echo off

setlocal enabledelayedexpansion
SET PARENT_DIR=%~dp0
SET PARENT_DIR=%PARENT_DIR:\=/%

if not defined SP_DIR (
    set SP_DIR=%CONDA_PREFIX%\Lib\site-packages
)
if not defined LIBRARY_PREFIX (
    set LIBRARY_PREFIX=%CONDA_PREFIX%\Library
)
if not defined LIBRARY_BIN (
    set LIBRARY_BIN=%LIBRARY_PREFIX%\bin
)

if not defined LIBRARY_LIB (
    set LIBRARY_LIB=%LIBRARY_PREFIX%\lib
)

REM Move code_aster and run_aster directories (including subdirectories)
echo Moving code_aster and run_aster directories...
move "%LIBRARY_PREFIX%\lib\aster\code_aster" "%SP_DIR%\code_aster"
move "%LIBRARY_PREFIX%\lib\aster\run_aster" "%SP_DIR%\run_aster"

REM Move all .pyd files to %SP_DIR%
echo Moving .pyd files...
for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.pyd") do move "%%f" "%SP_DIR%"

REM Move all dll/pdb files to %LIBRARY_BIN%
echo Moving .dll and .pdb files...
for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.dll") do move "%%f" "%LIBRARY_BIN%"
for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.pdb") do move "%%f" "%LIBRARY_BIN%"

REM Move all lib files to %LIBRARY_DIR%
echo Moving .lib files...
for %%f in ("%LIBRARY_PREFIX%\lib\aster\*.lib") do move "%%f" "%LIBRARY_LIB%"

REM if %CONDA_PREFIX%/Library/share/aster/tests not exists, create it
if not exist "%CONDA_PREFIX%\Library\share\aster\tests" (
    echo Creating %CONDA_PREFIX%/Library/share/aster/tests... from %PARENT_DIR%..\..\astest
    xcopy "%PARENT_DIR%..\..\astest\*" "%CONDA_PREFIX%\Library\share\aster\tests" /E /I /Y
)

echo Files moved successfully.

endlocal