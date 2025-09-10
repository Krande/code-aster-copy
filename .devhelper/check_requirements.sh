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

echo "code_aster requirements not found in 'build/'."
echo "If they are already installed, define ASTER_CONFIG environment variable" \
     "to the environment file."
printf "✨ Do you want to download the requirements ([y]/n)? "
read answer
answer=$(sed 's/^ *[nN] *$/no/' <<< "${answer}")
if [ "${answer}" = "no" ]; then
    exit 0
fi


_test curl || exit 1

is_intranet
okintranet=$?

mkdir -p ${prefix}/build
echo
if [ ${okintranet} -eq 0 ]; then
    echo "⏳ downloading requirements archive (it may take a few minutes)..."
    curl -fsSL https://minio.retd.edf.fr/codeaster/prereq/codeaster-prerequisites-${VERSION}-full-bin.sh \
        -o ${prefix}/build/archive.sh
else
    echo "❌ downloading requirements from internet is not yet supported, sorry"
    exit 1
fi
if [ ! -f ${prefix}/build/archive.sh ]; then
    echo "❌ download failed"
    exit 1
fi

cd ${prefix}/build
chmod 755 archive.sh
./archive.sh && rm -f archive.sh

printf "✅ code_aster requirements installed\n\n"

exit 0
