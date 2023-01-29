conda activate waf-ca-build && waf --prefix=%CONDA_PREFIX% --check-c-compiler=gcc --check-cxx-compiler=g++ configure install
rem conda activate waf-ca-build && waf --prefix=%CONDA_PREFIX% --without-hg --enable-metis --embed-metis --enable-mumps --embed-mumps --install-tests --disable-mfront --disable-petsc configure install --check-c-compiler=gcc --check-cxx-compiler=g++
rem ./waf build -j $CPU_COUNT
rem ./waf install