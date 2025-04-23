#!/bin/bash

# NB: 'profile.bat' should reflect the prerequisites paths.
# Paths are currently hard-coded in '{run_aster,run_ctest}.bat'...

if [ $# -ne 1 ]; then
    echo "usage: post_build_win.sh installdir"
    exit 1
fi

prefix="$1"

rm -rf ${prefix}/Python37
rm -rf ${prefix}/tools
rm -rf ${prefix}/med
rm -rf ${prefix}/medcoupling

cp -a /opt/public/win/Python37 ${prefix}/
cp -a /opt/public/win/tools ${prefix}/
cp -a /opt/public/win/med-4.1.1 ${prefix}/med
cp -a /opt/public/win/MEDCOUPLING_9_11_0 ${prefix}/medcoupling
