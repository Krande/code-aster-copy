#!/bin/bash -e

env

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
DEVTOOLS_URL=${ROOT_URL}/devtools.git
git clone ${DEVTOOLS_URL} devtools
(cd devtools ; git checkout main)

echo "+ downloading data..."
DATA_URL=${ROOT_URL}/data.git
git clone ${DATA_URL} data-src
(
    cd data-src
    branch=${CI_COMMIT_REF_NAME}
    echo "+ trying to fetch branch: ${branch}"
    git fetch origin ${branch}:${branch} > /dev/null 2>&1 || branch=v16
    echo "+ checking out branch: ${branch}"
    git checkout ${branch}
    git rev-parse --verify ${branch}
)
