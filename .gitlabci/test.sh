#!/bin/bash

jobs=${NPROC_MAX}
args=( "--clean" "--timefactor=4.0" "--jobs=${jobs}" "$@" )
# to be directly readable in the browser
export MESS_EXT="mess.txt"

if [ ! -z "${GITLAB_CI}" ]; then
    echo "+ fetching '${REFREV}' branch..."
    git branch -D ${REFREV} || true
    git fetch --depth=50 origin ${REFREV}
    git branch ${REFREV} FETCH_HEAD
fi
base=$(git merge-base ${REFREV} HEAD)

# check if it only changed testcases
changes=$(git diff --name-status ${base} | grep -v astest/)
if [ -z "${changes}" ]; then
    files=$(git diff --name-status ${base} | awk '{print $2}' | \
        grep astest/ | sed -e 's%astest/%%')
    if [ -z "${files}" ]; then
        printf "no changes detected?!\n"
    else
        ftmp=$(mktemp tmp.list.XXXXXXXX)
        for name in "${files}"; do
            grep -l "${name}" astest/*.export >> ${ftmp}
        done
        flist=$(mktemp tmp.list.XXXXXXXX)
        sort -u ${ftmp} | sed -e 's%astest/%%' -e 's%\.export%%' > ${flist}
        printf "\nchanged testcases:\n"
        cat ${flist}

        # only run these testcases
        args+=( "--testlist=${flist}" )
    fi
fi

# keep only outputs for failed tests, except for scheduled runs
if [ "${CI_PIPELINE_SOURCE}" != "schedule" ]; then
    args+=( "--only-failed-results" )
fi

printf "\nrun_ctest arguments: ${args}\n"

printf "\nrunning testcases #1... - $(date)\n"
./install/bin/run_ctest "${args[@]}"
iret=$?

if [ ${iret} -ne 0 ]; then
    printf "\nrunning testcases #2 (rerun-failed)... - $(date)\n"
    ./install/bin/run_ctest "${args[@]}" --rerun-failed
    iret=$?
fi

if [ ${iret} -ne 0 ]; then
    printf "\nrunning testcases #3 (rerun-failed)... - $(date)\n"
    ./install/bin/run_ctest "${args[@]}" --rerun-failed
    iret=$?
fi

# scheduled runs: archive results files
if [ "${CI_PIPELINE_SOURCE}" = "schedule" ]; then
    cd results
    tar czf mess_files.tar.gz *.${MESS_EXT}
    tar czf code_files.tar.gz *.code
    rm -f *.${MESS_EXT} *.code

    wget --no-verbose --no-check-certificate -O ./mc ${MINIO_URL}/codeaster/tools/mc
    chmod 755 ./mc
    ./mc --insecure alias set minio/ ${MINIO_URL} ${MINIO_LOGIN} ${MINIO_PASSWD}
    dest=minio/codeaster/devops/ci-debian11/results/${REFREV}/verification
    ./mc --insecure cp run_testcases.xml ${dest}/
    ./mc --insecure cp mess_files.tar.gz ${dest}/
    ./mc --insecure cp code_files.tar.gz ${dest}/

    cd ..
fi

exit ${iret}
