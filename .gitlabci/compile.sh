#!/bin/bash -e

echo "+ compiling..."

opts=( "--prefix=./install" "--without-repo" )
[ -f data-src/README ] && opts+=( "--with-data=data-src" )

./configure "${opts[@]}"

jobs=$(( ${NPROC_MAX} / 2 ))
make install -j ${jobs}
