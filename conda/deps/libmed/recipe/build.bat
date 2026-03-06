@ECHO ON

setlocal

mkdir build
cd build

:: Needed by IFX
set "LIB=%BUILD_PREFIX%\Library\lib;%LIB%"
set "INCLUDE=%BUILD_PREFIX%\opt\compiler\include\intel64;%INCLUDE%"

set FCFLAGS=/fpp /nologo %FCFLAGS%

cmake -G "Ninja" ^
  %CMAKE_ARGS% ^
  -D HDF5_BUILD_FORTRAN:BOOL=ON ^
  -D CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS:BOOL=ON ^
  -D Python_FIND_STRATEGY:STRING=LOCATION ^
  -D Python_FIND_REGISTRY:STRING=NEVER ^
  -D Python3_ROOT_DIR:FILEPATH="%PREFIX%" ^
  -D HDF5_ROOT_DIR:FILEPATH="%LIBRARY_PREFIX%" ^
  -D CMAKE_Fortran_FLAGS:STRING="%FCFLAGS%" ^
  -D MEDFILE_INSTALL_DOC=OFF ^
  -D MEDFILE_BUILD_PYTHON=ON ^
  -D MEDFILE_BUILD_TESTS=OFF ^
  -D MEDFILE_BUILD_SHARED_LIBS=ON ^
  -D MEDFILE_BUILD_STATIC_LIBS=OFF ^
  -D MEDFILE_USE_UNICODE=OFF ^
  -D MED_MEDINT_TYPE=int ^
  ..

if errorlevel 1 exit 1

ninja
if errorlevel 1 exit 1

mkdir %SP_DIR%\med
if errorlevel 1 exit 1
ninja install
if errorlevel 1 exit 1

:: Regenerate import libraries with lowercase+underscore aliases (e.g. mfacre_ = MFACRE)
:: so that code_aster compiled with /names:lowercase /assume:underscore can link against libmed
python %RECIPE_DIR%\fix_exports.py %LIBRARY_BIN% %LIBRARY_LIB%
if errorlevel 1 exit 1

copy %LIBRARY_BIN%\mdump4.exe %LIBRARY_BIN%\mdump.exe
if errorlevel 1 exit 1
copy %LIBRARY_BIN%\xmdump4 %LIBRARY_BIN%\xmdump
if errorlevel 1 exit 1

endlocal