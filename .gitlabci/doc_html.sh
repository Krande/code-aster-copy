#!/bin/bash -e

printf "\nGenerate html documentation...\n"
make doc

getstatus="git status -uno --porcelain"
if [ $(${getstatus} | wc -l) != 0 ]; then
    printf "\nChanges must be committed:\n"
    ${getstatus}
    exit 1
fi
