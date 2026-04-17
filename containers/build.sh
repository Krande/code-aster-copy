#!/bin/bash -ex
# because it may use dash...

mkdir -p ${WRKDIR}
cd ${WRKDIR}

echo "+ cloning..."
git clone --depth=1 --branch=${BRANCH} https://gitlab.pleiade.edf.fr/codeaster/src.git
cd ${WRKDIR}/src

echo "+ creating pkginfo..."
echo "${PKGINFO}" > code_aster/pkginfo.py
cat code_aster/pkginfo.py

echo "+ running configure..."
source /opt/public/${DEVTOOLS_COMPUTER_ID}_mpi.sh
./configure --prefix=/opt/codeaster

echo "+ compiling..."
jobs=$(( $(nproc) - 2 ))
make install -j ${jobs}

echo "+ cleaning..."
cd /
rm -rf ${WRKDIR}
