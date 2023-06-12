#!/bin/bash -e

echo "+ compiling..."
export CONFIG_PARAMETERS_mpiexec="mpiexec -n {mpi_nbcpu} --allow-run-as-root --tag-output {program}"
export DEVTOOLS_COMPUTER_ID=none

# only add mpi4py
. env.d/version.sh
export PREREQ_PATH=/opt/public/${VERSION}/gcc8-openblas-ompi3
export PYPATH_MPI4PY="${PREREQ_PATH}/mpi4py-3.1.3/lib/python3.7/site-packages"
export PYTHONPATH="${PYPATH_MPI4PY}:${PYTHONPATH}"

# debug build
export BUILD=debug
./configure --prefix=./mini --without-repo --no-enable-all

jobs=$(( $(nproc) / 2 ))
make install -j ${jobs}
