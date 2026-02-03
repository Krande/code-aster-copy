#!/bin/bash -e

set_prefix() {
    local this=$(readlink -n -f "$1")
    BASE=$(dirname "${this}")
}

set_prefix "${0}"

MINIO_DIR=minio/codeaster/devops/perf
host=${OSNAME:-$(uname -n)}

_help() {
    echo usage: $(basename $0) results-directory upload
}

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    _help
    exit 1
fi

resdir=$1
upload=$2
wrkdir=build/perf
mkdir -p ${wrkdir}

printf "\nperf.sh - start time - $(date)\n"

if [ ! -z "${MINIO_URL}" ]; then
    printf "\nsetting up minIO alias - $(date)\n"
    mc --insecure alias set minio/ ${MINIO_URL} ${MINIO_LOGIN} ${MINIO_PASSWD}
fi

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
    printf "code files are supposed to be here (ci job)\n"
fi

today=$(date +%Y-%m-%d)
printf "\nextracting data from '.code' files...\n"
python3 .devhelper/check_perf.py --extract -o ${wrkdir}/${today}.csv ${resdir}

printf "\nsearching for previous data... - $(date)\n"
csv=( $(mc ls --json ${MINIO_DIR}/${host} | jq -r '"\(.key) \(.lastModified)"') )
idx=0
size=$(( ${#csv[@]} / 2 ))
prev=( $(
    (
    for ((i=0;i<${size};i++)); do
        in=$(( 2 * $i ))
        idate=$(( 2 * $i + 1 ))
        filename=${csv[$in]}
        date=${csv[$idate]:0:10}
        if [ "${filename##*.}" != "csv" ]; then
            continue
        fi
        echo ${date} ${filename}
    done
    ) | sort | tail -1
) )

printf "\ndownload previous data into ${wrkdir} - $(date)\n"
# 'prev' contains: last-date last-filename
mc --insecure cp ${MINIO_DIR}/${host}/${prev[1]} ${wrkdir}/${prev[0]}.csv

last=${wrkdir}/last_changes_${host}.txt
if [ ${#prev[@]} -eq 0 ]; then
    printf "\nWARNING no previous data for comparison.\n\n"
else
    printf "\nchecking changes since previous run...\n"
    python3 .devhelper/check_perf.py --compare ${wrkdir}/${today}.csv ${wrkdir}/${prev[0]}.csv > ${last}

    printf "\n----------------------------------\n"
    cat ${last}
    printf "\n----------------------------------\n"
fi

printf "\nupload analysis onto minio (upload=${upload}) - $(date)\n"
if [ "${upload}" = "upload" ]; then
    mc --insecure cp ${wrkdir}/${today}.csv ${MINIO_DIR}/${host}/
    [ -f ${last} ] && mc --insecure cp ${last} ${MINIO_DIR}/
else
    echo "INFO: do not upload results"
fi

printf "\nend time - $(date)\n\n"
exit 0
