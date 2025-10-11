# This file set the environment for code_aster.
# Configuration for cronos mpi
. $(readlink -n -f $(dirname ${BASH_SOURCE:-${(%):-%x}}))/version.sh
. /software/restricted/simumeca/aster/prerequisites/${VERSION}/gcc8-mkl-ompi4/cronos_mpi.sh
