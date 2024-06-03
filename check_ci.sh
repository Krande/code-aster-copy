#!/bin/bash

# shell script to simulate (some steps of) the CI procedure

export SUBMIT_MODE=1

args=( "${@}" )
if [ $# -eq 0 ]; then
    args=( "--all" )
fi

installdir=${PREFIX_ROOT:-../install}
ln -sf ../devtools .
ln -sf ${installdir}/mpi install
trap "rm -f devtools install" EXIT

.gitlabci/debug-ci.sh "${args[@]}"
