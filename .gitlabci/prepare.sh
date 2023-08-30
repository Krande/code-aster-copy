#!/bin/bash -e

env
# CI_PROJECT_URL=https://gitlab.pleiade.edf.fr/codeaster/lab/experiment/src
root=$(dirname ${CI_PROJECT_URL})
if grep -q experiment <<< "${CI_PROJECT_URL}"; then
    root=$(dirname $(dirname ${root}))
fi

echo "+ downloading the runner image..."
source env.d/version.sh
URL_SIF=${MINIO_URL}/codeaster/sif/ci/codeaster-prerequisites-${VERSION}-debian-10.sif
if [ ! -z ${SIF} ]; then
    if [ -z "${DEBUG_CI}" ] || [ ! -f ${ORIG_HOME}/containers/$(basename ${URL_SIF}) ]; then
        wget --no-check-certificate -O ${SIF} ${URL_SIF}
    else
        cp ${ORIG_HOME}/containers/$(basename ${URL_SIF}) ${SIF}
    fi
fi

echo "+ downloading devtools..."
DEVTOOLS_URL=${root}/devtools
grep -q https: <<< "${DEVTOOLS_URL}" && DEVTOOLS_URL=${DEVTOOLS_URL}.git
git clone ${DEVTOOLS_URL} devtools
(cd devtools ; git checkout main)

echo "+ downloading data..."
DATA_URL=${root}/data
grep -q https: <<< "${DATA_URL}" && DATA_URL=${DATA_URL}.git
git clone ${DATA_URL} data-src
(
    cd data-src
    branch=${CI_COMMIT_BRANCH}
    git rev-parse --verify ${branch} > /dev/null 2>&1 || branch=v15
    git checkout ${branch}
)
