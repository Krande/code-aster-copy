#!/bin/bash -e

env

echo "+ checking for requirements..."
source env.d/version.sh
if [ "${VERSION}" != "${PREREQ_VERSION}" ]; then
    echo "Docker image (${PREREQ_VERSION}) and prerequisites version (${VERSION}) are inconsistent"
    exit 1
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
    git fetch origin ${branch}:${branch} > /dev/null 2>&1 || branch=${REFREV}
    echo "+ checking out branch: ${branch}"
    git checkout ${branch}
    git rev-parse --verify ${branch}
)
