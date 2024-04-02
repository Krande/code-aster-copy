#!/bin/bash
export -e
export CLICOLOR_FORCE=1

export PREFIX="${CONDA_PREFIX}"
echo "PREFIX=${PREFIX}"

#export LD_LIBRARY_PATH=${PREFIX}/lib:$LD_LIBRARY_PATH
#export LDFLAGS=${PREFIX}/lib:${LDFLAGS}

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
export INCLUDES_MUMPS=${PREFIX}/include

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
     --med-libs=medC \
     --mumps-libs="dmumps_seq zmumps_seq smumps_seq cmumps_seq mumps_common_seq pord_seq" \
     --install-tests \
     --maths-libs=auto \
     --disable-mpi \
     --without-hg \
     --without-repo \
     --conda-build \
     configure

./waf_std install_debug -v
