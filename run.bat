conda activate waf-ca-build &&^
waf distclean && ^
waf configure^
    --prefix=%CONDA_PREFIX% ^
    --includedir=%CONDA_PREFIX%\include ^
    --libdir=%CONDA_PREFIX%\lib ^
    -kk ^
    --disable-openmp ^
    --shared-aster ^
    --disable-mfront ^
    --disable-petsc ^
    --check-c-compiler=gcc^
    --check-cxx-compiler=g++^
    --without-hg && ^
waf install
rem ./waf build -j $CPU_COUNT
rem ./waf install