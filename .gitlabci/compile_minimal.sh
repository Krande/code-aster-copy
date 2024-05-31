#!/bin/bash -e

echo "+ compiling..."
export DEVTOOLS_COMPUTER_ID=none

# only add mpi4py
. env.d/version.sh
export PREREQ_PATH=/opt/public/${VERSION}/gcc-openblas-ompi
export PYPATH_MPI4PY="$(find ${PREREQ_PATH}/mpi4py-*/lib/python* -name site-packages)"
export PYTHONPATH="${PYPATH_MPI4PY}:${PYTHONPATH}"

jobs=$(( ${NPROC_MAX} / 2 ))
export BUILD=debug

# mpi build
./configure --prefix=./mini --without-repo --no-enable-all
make install -j ${jobs}
make distclean

# sequential build
./configure --prefix=./mini --without-repo --no-enable-all --disable-mpi
make install -j ${jobs}
