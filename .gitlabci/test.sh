#!/bin/bash

jobs=$(( $(nproc) * 7 / 8 ))
args=( "--clean" "--timefactor=4.0" "--jobs=${jobs}" "$@" )
# to be directly readable in the browser
export MESS_EXT="mess.txt"

refrev=v16

if [ ! -z "${GITLAB_CI}" ]; then
    echo "+ fetching '${refrev}' branch..."
    git branch -D ${refrev} || true
    git fetch --depth=50 origin ${refrev}
    git branch ${refrev} FETCH_HEAD
fi
base=$(git merge-base ${refrev} HEAD)

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

exit ${iret}
