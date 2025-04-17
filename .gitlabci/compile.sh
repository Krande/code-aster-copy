#!/bin/bash -e

echo "+ compiling..."

opts=( "--prefix=./install" "--without-repo" )
[ -f data-src/README ] && opts+=( "--with-data=data-src" )

if [ "${BUILDTYPE}" = "nightly-coverage" ]; then
    if [ "${ASTER_BUILD}" != "debug" ]; then
        echo "ERROR: ASTER_BUILD must be set as 'debug'"
        exit 4
    fi
    opts+=( "--coverage" )
fi

jobs=$(( ${NPROC_MAX} / 2 ))

if [ "${OSNAME}" != "win" ]; then
    ./configure "${opts[@]}"
    make install -j ${jobs}
else
    ./waf_std configure --mingw-cross-compilation "${opts[@]}"
    ./waf_std install -j ${jobs}

    cp -a /opt/public/win/Python37 .
    cp -a /opt/public/win/tools outils
    cp -a /opt/public/win/MEDCOUPLING_9_11_0 ./install/medcoupling
fi
