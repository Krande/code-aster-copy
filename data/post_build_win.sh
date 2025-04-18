#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: post_build_win.sh installdir"
    exit 1
fi

prefix="$1"
root="$prefix/.."

mkdir -p ${root}/outils
rm -rf ${root}/Python37
rm -rf ${prefix}/medcoupling

cp -a /opt/public/win/Python37 ${root}/
cp -a /opt/public/win/tools/* ${root}/outils/
cp -a /opt/public/win/MEDCOUPLING_9_11_0 ${prefix}/medcoupling
