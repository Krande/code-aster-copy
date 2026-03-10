@echo on

set src=%cd%

cd work
copy %RECIPE_DIR%\CMakeLists.txt %src%\CMakeLists.txt

:: Use 64-bit MUMPS_INT to match code_aster's /integer-size:64
copy %src%\src\mumps_int_def64_h.in %src%\include\mumps_int_def.h

:: Set up Intel Fortran compiler paths from conda build prefix
set "PATH=%BUILD_PREFIX%\Library\bin\compiler;%BUILD_PREFIX%\Library\bin;%BUILD_PREFIX%\Scripts;%PATH%"
set "LIB=%BUILD_PREFIX%\Library\lib;%LIB%"
set "INCLUDE=%BUILD_PREFIX%\opt\compiler\include\intel64;%BUILD_PREFIX%\Library\include;%INCLUDE%"
set FC=ifx

mkdir build
cd build

:: Configure: ILP64 mode with debug symbols
cmake -G "Ninja" ^
      -DCMAKE_Fortran_COMPILER=ifx ^
      -DCMAKE_PREFIX_PATH=%LIBRARY_PREFIX% ^
      -DCMAKE_INSTALL_PREFIX:PATH=%LIBRARY_PREFIX% ^
      -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo ^
      -DUSE_ILP64:BOOL=ON ^
      ..
if errorlevel 1 exit 1
cmake --build . --config RelWithDebInfo --target install
if errorlevel 1 exit 1

:: Verify simpletests
%src%\build\c_example < %src%\examples\input_simpletest_real
if errorlevel 1 exit 1
%src%\build\ssimpletest < %src%\examples\input_simpletest_real
if errorlevel 1 exit 1
%src%\build\dsimpletest < %src%\examples\input_simpletest_real
if errorlevel 1 exit 1
%src%\build\csimpletest < %src%\examples\input_simpletest_cmplx
if errorlevel 1 exit 1
%src%\build\zsimpletest < %src%\examples\input_simpletest_cmplx
if errorlevel 1 exit 1
