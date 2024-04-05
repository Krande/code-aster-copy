#!/bin/bash
set -e

which python
# set unbuffered output
export PYTHONUNBUFFERED=1

if [ -z "${CONDA_PREFIX}" ]; then
  # if CONDA_PREFIX is not defined, we are not in a conda environment and will have to activate it using the env var
  # CONDA_ENV_DIR and activate it. The CONDA_ENV_DIR is defined in the .env file
  for line in $(cat .env); do export $line; done
  if [ -z "${CONDA_ENV_DIR}" ]; then
      echo "CONDA_ENV_DIR is not defined"
      exit 1
  fi
  source ${CONDA_ENV_DIR}/bin/activate codeaster-deps
else
  echo "CONDA_PREFIX is defined"
  export CLICOLOR_FORCE=1
fi

export PREFIX="${CONDA_PREFIX}"
echo "PREFIX=${PREFIX}"

export LD_LIBRARY_PATH="${PREFIX}/lib ${LD_LIBRARY_PATH}"

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

export DEFINES="H5_BUILT_AS_DYNAMIC_LIB H5_USE_110_API ${DEFINES}"

which gcc
which g++
which gfortran

python conda/update_version.py

./waf distclean

./waf_std \
     --python=$PYTHON \
     --prefix="${PREFIX}" \
     --libdir="${PREFIX}/lib" \
     --enable-hdf5 \
     --enable-med \
     --med-libs="med medC medfwrap medimport" \
     --mumps-libs="dmumps_seq zmumps_seq smumps_seq cmumps_seq mumps_common_seq pord_seq" \
     --install-tests \
     --maths-libs="cblas lapack" \
     --disable-mpi \
     --without-hg \
     --without-repo \
     configure

./waf_std install_debug -v
