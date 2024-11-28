@echo off

setlocal enabledelayedexpansion


if defined OUTPUT_DIR (
    echo "OUTPUT_DIR: %OUTPUT_DIR%"
) else (
    set OUTPUT_DIR=temp\\seq
)

run_ctest --resutest=%OUTPUT_DIR% -L submit -L sequential -LE need_data --timefactor=10.0 --only-failed-results

if errorlevel 1 exit 1


