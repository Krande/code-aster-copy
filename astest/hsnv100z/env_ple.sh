#!/bin/bash

if [ ! -z "${SINGULARITY_NAME}" ]; then
    # In scibian container:
    # . /opt/public/debian-11_mpi.sh
    export CA_INSTALL=${HOME}/dev/codeaster/install/mpi
    export CS_INSTALL=${HOME}/dev/salome_cfd/cs_install/old
    PYV=3.9
else
    # On cronos:
    # . ${HOME}/dev/codeaster-git/src/env.d/cronos_mpi.sh
    export CA_INSTALL=/software/restricted/simumeca/aster/install/unstable/mpi
    export CS_INSTALL=/software/rd/saturne/code_saturne/nightly/arch/cronos_ast_coupling_dbg
    PYV=3.6
fi

# add pyple path
LIBPLE_DIR=${CS_INSTALL}

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LIBPLE_DIR}/lib
export PYTHONPATH=${PYTHONPATH}:${LIBPLE_DIR}/lib/python${PYV}/site-packages

# for code_aster runner
export PATH=${CA_INSTALL}/bin:${PATH}
