#!/bin/bash

set_prefix() {
    local this=$(readlink -n -f "$1")
    BASE=$(dirname "${this}")
}

set_prefix "${0}"

MINIO_DIR=minio/codeaster/devops/coverage

_help() {
    echo usage: $(basename $0) src-directory install-directory results-directory upload
}

if [ $# -lt 3 ] || [ $# -gt 4 ]; then
    _help
    exit 1
fi

srcdir=$1
instdir=$2
resdir=$3
upload=$4
wrkdir=build/coverage

printf "\ncoverage.sh - start time - $(date)\n"

if [ ! -z "${MINIO_URL}" ]; then
    printf "\nsetting up minIO alias - $(date)\n"
    mc --insecure alias set minio/ ${MINIO_URL} ${MINIO_LOGIN} ${MINIO_PASSWD}
fi

printf "\ndownload previous data into ${wrkdir} - $(date)\n"
mkdir -p ${wrkdir}
mc --insecure cp ${MINIO_DIR}/last.tested ${wrkdir}/

printf "\nchecking testcases results directory ${resdir} - $(date)\n"
mkdir -p ${resdir}

if [ -f ${resdir}/code_files.tar.gz ]; then
    (
        cd ${resdir}
        printf "extracting code files...\n"
        tar xzf code_files.tar.gz
        printf "files extracted.\n"
    )
else
    printf "code files supposed to be here (merge-request job)\n"
fi

printf "\ncoverage analysis - $(date)\n"
source ${instdir}/share/aster/profile.sh
python3 ${BASE}/coverage_main.py --verbose \
    --srcdir=${srcdir} --installdir=${instdir} --resdir=${resdir} \
    --wrkdir=${wrkdir} --previous=last \
    --save --savetxt
iret=$?
# --limit=180

printf "\nupload analysis onto minio (upload=${upload}) - $(date)\n"
if [ "${upload}" = "upload" ]; then
    today=$(date +%Y-%m-%d)
    mc --insecure cp ${wrkdir}/${today}.* ${MINIO_DIR}/
    mc --insecure cp ${wrkdir}/last.* ${MINIO_DIR}/
else
    echo "INFO: do not upload results"
fi

printf "\nend time - $(date)\n\n"
exit ${iret}
