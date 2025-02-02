@echo off
echo "Activating Intel environment..."
call "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env\vars.bat" -arch intel64
if errorlevel 1 (
    echo "Failed to activate Intel environment."
    exit /b 1
)
where ifx.exe
if errorlevel 1 exit /b 1
