# This file set the environment for code_aster.
# Configuration for selena mpi
. $(readlink -n -f $(dirname ${BASH_SOURCE:-${(%):-%x}}))/version.sh

. /software/shared/simumeca/aster/prerequisites/${VERSION}/gcc-mkl-ompi509/selena_mpi.sh
