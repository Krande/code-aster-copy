#!/bin/bash

set -x
variant="$1"

echo "+ $(date) - building code_aster Docker image on ${variant}..."
baseimage="${DOCKER_NEXUS_URL}/codeaster-prerequisites:${PREREQ_VERSION}-${variant}"
dest="${DOCKER_NEXUS_URL}/codeaster-main:${PREREQ_VERSION}-${variant}"

echo "${DOCKER_NEXUS_PASSWD}" | docker login ${DOCKER_NEXUS_URL} -u "${DOCKER_NEXUS_USER}" --password-stdin

cd .gitlabci
docker build \
    --build-arg baseimage=${baseimage} \
    --build-arg osname=${variant} \
    -t ${dest} .

docker image push ${dest}
