#!/bin/bash
set -e

export CLICOLOR_FORCE=1

which python

export PREFIX="${CONDA_PREFIX}"
echo "PREFIX=${PREFIX}"
# If prefix is not set, activate CONDA environment
#export LD_LIBRARY_PATH="${PREFIX}/lib ${LD_LIBRARY_PATH}"
#export LDFLAGS=${PREFIX}/lib:${LDFLAGS}
echo "LDFLAGS"
#export LIBPATH="$PREFIX/lib $LIBPATH"
#export LDFLAGS="-Wl,--no-as-needed -L$PREFIX/lib -lm -lpthread -ldl -lz -lgomp ${LDFLAGS}"

export CC=gcc
export CXX=g++
export FC=gfortran

export LIBPATH_HDF5=${PREFIX}/lib
export INCLUDES_HDF5=${PREFIX}/include

export LIBPATH_MED=${PREFIX}/lib
export INCLUDES_MED=${PREFIX}/include

export LIBPATH_METIS=${PREFIX}/lib
export INCLUDES_METIS=${PREFIX}/include

export LIBPATH_MUMPS=${PREFIX}/lib
export INCLUDES_MUMPS="${PREFIX}/include ${PREFIX}/include/mumps_seq"

export LIBPATH_SCOTCH=${PREFIX}/lib
export INCLUDES_SCOTCH=${PREFIX}/include

export TFELHOME=${PREFIX}

export LIBPATH_MGIS=${PREFIX}/bin
export INCLUDES_MGIS=${PREFIX}/include

export FCFLAGS="-fallow-argument-mismatch ${FCFLAGS}"

which gcc
which g++
which gfortran

./waf distclean

./waf_std \
     --python=$PYTHON \
     --prefix="${PREFIX}" \
     --med-libs="med medC medfwrap medimport" \
     --mumps-libs="dmumps_seq zmumps_seq smumps_seq cmumps_seq mumps_common_seq pord_seq" \
     --install-tests \
     --maths-libs=auto \
     --disable-mpi \
     --without-hg \
     --without-repo \
     configure

./waf_std install_debug -v
