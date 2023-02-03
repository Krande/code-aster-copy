#!/bin/bash -e

echo "+ compiling..."
export CONFIG_PARAMETERS_mpiexec="mpiexec -n {mpi_nbcpu} --allow-run-as-root --tag-output {program}"

# --with-data=data-src
./configure --prefix=./install --without-repo
make install
