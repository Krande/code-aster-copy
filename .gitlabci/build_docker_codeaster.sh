#!/busybox/sh
# Warning: will be executed by /busybox/sh in a kaniko container

variant="$1"

echo "+ $(date) - building code_aster Docker image on ${variant}..."
baseimage="${DOCKER_NEXUS_URL}/codeaster-prerequisites:${PREREQ_VERSION}-${variant}"
dest="${DOCKER_NEXUS_URL}/codeaster-main:${PREREQ_VERSION}-${variant}"

# instanciate template
mkdir -p artf
sed -e "s@%baseimage%@${baseimage}@g" .gitlabci/codeaster-main.dockerfile.tmpl \
    > artf/codeaster-main-${variant}.dockerfile
cat artf/codeaster-main-${variant}.dockerfile

encoded=$(printf "%s:%s" "${DOCKER_NEXUS_USER}" "${DOCKER_NEXUS_PASSWD}" | base64 | tr -d '\n')
cat << eof > /kaniko/.docker/config.json
{
    "auths": {
        "${DOCKER_NEXUS_URL}": {
            "auth": "${encoded}"
        }
    }
}
eof

# copying certificate into /kaniko/ssl/certs/ca-certificates.crt does not work
/kaniko/executor \
      --context "${CI_PROJECT_DIR}" \
      --build-arg GITLAB_PREREQ_TOKEN=${GITLAB_PREREQ_TOKEN} \
      --dockerfile "${CI_PROJECT_DIR}/artf/codeaster-main-${variant}.dockerfile" \
      --skip-tls-verify \
      --destination "${dest}"
