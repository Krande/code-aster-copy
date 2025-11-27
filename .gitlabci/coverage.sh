#!/bin/bash

set_prefix() {
    local this=$(readlink -n -f "$1")
    BASE=$(dirname "${this}")
}

set_prefix "${0}"

MINIO_DIR=minio/codeaster/devops/coverage

_help() {
    echo usage: $(basename $0) src-directory install-directory results-directory
}

if [ $# -ne 3 ]; then
    _help
    exit 1
fi

srcdir=$1
instdir=$2
resdir=$3
wrkdir=build/coverage
mkdir -p ${wrkdir}

printf "\nstart time : $(date)\n"
printf "\ndownload previous data into ${wrkdir} : $(date)\n"
mc --insecure cp ${MINIO_DIR}/last.tested ${wrkdir}/

printf "\ncoverage analysis : $(date)\n"
source ${instdir}/share/aster/profile.sh
python3 ${BASE}/coverage_main.py --verbose \
    --srcdir=${srcdir} --installdir=${instdir} --resdir=${resdir} \
    --wrkdir=${wrkdir} --previous=last --limit=180 \
    --save --savetxt

printf "\nupload analysis onto minio : $(date)\n"
if [ -z "${CI_COMMIT_REF_NAME}" ]; then
    echo "WARNING: results are only uploaded in CI"
else
    today=$(date +%Y-%m-%d)
    mc --insecure ${wrkdir}/${today}.* ${MINIO_DIR}/
    mc --insecure ${wrkdir}/last.* ${MINIO_DIR}/
fi

printf "end time : $(date)\n\n"
