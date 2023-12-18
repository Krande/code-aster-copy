#!/bin/bash

# shell script to simulate (some steps of) the CI procedure

export SUBMIT_MODE=1

ln -sf ../devtools .
ln -sf ~/.install/debian-10/mpi install
trap "rm -f devtools install" EXIT

args=( "${@}" )
if [ $# -eq 0 ]; then
    args=( "--all" )
fi

.gitlabci/debug-ci.sh "${args[@]}"
