#!/bin/bash

prefix=$(dirname $(dirname $(readlink -n -f ${0})))

_test() {
    printf "\nchecking for ${1}...\n"
    ${1} --version 2> /dev/null
    if [ $? -ne 0 ]; then
        echo "ERROR: please install ${1}."
        return 1
    fi
    return 0
}

is_intranet() {
    printf "checking for EDF network... "
    curl -m 2 https://gitlab.pleiade.edf.fr > /dev/null 2>&1
    okintranet=$?
    test ${okintranet} -eq 0 && echo ok || echo no
    return ${okintranet}
}

echo
echo "code_aster requirements not found in 'build/'."
echo "If they are already installed, define ASTER_CONFIG environment variable" \
     "to the environment file."
printf "✨ Do you want to download the requirements ([y]/n, c to cancel)? "
read answer
answer=$(tr '[:upper:]' '[:lower:]' <<< "${answer}")
[ "${answer}" = "no" ] && exit 0
[ "${answer}" = "c" ] && echo "exiting..." && exit 1

ASTER_REQS_PACKAGE=${ASTER_REQS_PACKAGE:-"gcc-ompi"}
if [ "${ASTER_REQS_PACKAGE}" = "ask" ]; then
    echo
    echo "Select the package you want to download:"
    echo "  1. full (embedding gcc & openmpi, recommended)"
    echo "  2. with gcc, without openmpi (you must have the same version of openmpi on the host)"
    echo "  3. without gcc, without openmpi (you must have the same version of openmpi on the host)"
    printf "✨ your choice ([1]/2/3, 0 to cancel)? "
    read answer
    [ "${answer}" = "2" ] && ASTER_REQS_PACKAGE="gcc-noompi"
    [ "${answer}" = "3" ] && ASTER_REQS_PACKAGE="nogcc-noompi"
    [ "${answer}" = "0" ] && echo "exiting..." && exit 1
fi
arch="codeaster-prerequisites-${VERSION}-${ASTER_REQS_PACKAGE}.sh"

_test curl || exit 1

is_intranet
okintranet=$?

mkdir -p ${prefix}/build
echo
if [ ${okintranet} -eq 0 ]; then
    echo "⏳ downloading requirements archive ${arch} (it may take a few minutes)..."
    curl -fsSL https://minio.retd.edf.fr/codeaster/prereq/${arch} \
        -o ${prefix}/build/${arch}
else
    echo "❌ downloading requirements from internet is not yet supported, sorry"
    exit 1
fi
if [ ! -f ${prefix}/build/${arch} ]; then
    echo "❌ download failed"
    exit 1
fi

cd ${prefix}/build
chmod 755 ${arch}
./${arch} && rm -f ${arch}

printf "✅ code_aster requirements installed\n\n"

exit 0
