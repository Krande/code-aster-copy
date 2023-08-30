#!/bin/bash

jobs=$(( $(nproc) * 7 / 8 ))
args=( "--clean" "--only-failed-results" "--timefactor=4.0" "--jobs=${jobs}" "$@" )

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
