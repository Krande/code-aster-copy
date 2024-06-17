:: created by waf using data/wscript

:: IF ASTER_ROOT is not defined and CONDA_PREFIX is defined, use CONDA_PREFIX

if not defined ASTER_ROOT if defined CONDA_PREFIX set ASTER_ROOT=%CONDA_PREFIX%

set PYTHONPATH=%ASTER_ROOT%\Library\lib\aster;%ASTER_ROOT%\lib\site-packages;%PYTHONPATH%

set ASTER_DATADIR=%ASTER_ROOT%\Library\share\aster
set ASTER_LIBDIR=%ASTER_ROOT%\Library\lib\aster
set ASTER_LOCALEDIR=%ASTER_ROOT%\Library\share\locale\aster
set ASTER_ELEMENTSDIR=%ASTER_ROOT%\Library\lib\aster