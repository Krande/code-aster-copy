#!/bin/bash -e

env
echo "+ downloading the runner image..."
source env.d/version.sh
URL_SIF=${MINIO_URL}/codeaster/sif/ci/codeaster-prerequisites-${VERSION}-debian-10.sif
if [ ! -z ${SIF} ]; then
    if [ -z "${DEBUG_CI}" ] || [ ! -d ${HOME}/containers/$(basename ${URL_SIF}) ]; then
        wget --no-check-certificate -O ${SIF} ${URL_SIF}
    else
        cp ${HOME}/containers/$(basename ${URL_SIF}) ${SIF}
    fi
fi

# CI_REPOSITORY_URL=https://gitlab.pleiade.edf.fr/codeaster/lab/experiment/src
root=$(dirname ${CI_REPOSITORY_URL})
if grep -q experiment <<< "${CI_REPOSITORY_URL}"; then
    root=$(dirname $(dirname ${root}))
fi

echo "+ downloading devtools..."
DEVTOOLS_URL=${root}/devtools
grep -q https: <<< "${DEVTOOLS_URL}" && DEVTOOLS_URL=${DEVTOOLS_URL}.git
git clone ${DEVTOOLS_URL} devtools
(cd devtools ; git checkout use-git)

echo "+ downloading data..."
DATA_URL=${root}/data
grep -q https: <<< "${DATA_URL}" && DATA_URL=${DATA_URL}.git
git clone ${DATA_URL} data-src
(
    cd data-src
    branch=${CI_COMMIT_BRANCH}
    git rev-parse --verify ${branch} > /dev/null 2>&1 || branch=main
    git checkout ${branch}
)
