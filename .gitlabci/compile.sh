#!/bin/bash -e

echo "+ compiling..."

opts=( "--prefix=./install" "--without-repo" )
[ -f data-src/README ] && opts+=( "--with-data=data-src" )

if [ "${BUILDTYPE}" = "nightly-coverage" ]; then
    if [ "${BUILD}" != "debug" ]; then
        echo "ERROR: BUILD must be set as 'debug'"
        exit 4
    fi
    opts+=( "--coverage" )
fi

./configure "${opts[@]}"

jobs=$(( ${NPROC_MAX} / 2 ))
make install -j ${jobs}
