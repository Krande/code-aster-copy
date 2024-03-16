IF "%CONDA_BUILD%" == "" (
    call build_env_setup.bat
)
REM ===== end generated header =====
@echo off

REM Install for standard sequential
waf configure ^
  --prefix=%LIBRARY_PREFIX% ^
  --libdir=%LIBRARY_PREFIX%\lib ^
  --pythondir=%LIBRARY_PREFIX%\lib ^
  --install-tests ^
  --embed-metis ^
  --without-hg ^
  
waf install

