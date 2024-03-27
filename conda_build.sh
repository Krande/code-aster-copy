#!/bin/bash
set -e
export CLICOLOR_FORCE=1

export PREFIX="${CONDA_PREFIX}"
echo "PREFIX=${PREFIX}"

export LD_LIBRARY_PATH=${PREFIX}/lib:$LD_LIBRARY_PATH

export CC=gcc
export CXX=g++
export FC=gfortran


export FCFLAGS="-fallow-argument-mismatch ${FCFLAGS}"

which gcc
which g++
which gfortran

./waf_std \
     --python=$PYTHON \
     --prefix="${PREFIX}" \
     --libdir="${PREFIX}/lib" \
     --install-tests \
     --maths-libs=auto \
     --disable-mpi \
     --disable-petsc \
     --without-hg \
     --without-repo \
     configure

./waf_std install_debug
