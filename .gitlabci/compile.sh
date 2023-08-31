#!/bin/bash -e

echo "+ compiling..."
export CONFIG_PARAMETERS_mpiexec="mpiexec -n {mpi_nbcpu} --allow-run-as-root --tag-output {program}"

opts=( "--prefix=./install" "--without-repo" )
[ -f data-src/README ] && opts+=( "--with-data=data-src" )

. /opt/public/*_mpi.sh
./configure "${opts[@]}"

jobs=$(( $(nproc) / 2 ))
make install -j ${jobs}
