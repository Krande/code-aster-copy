# This file set the environment for code_aster.
# Configuration for selena mpi
. $(readlink -n -f $(dirname ${BASH_SOURCE:-${(%):-%x}}))/version.sh

# transitional
[ "${VERSION}" = "20240327" ] && VERSION="${VERSION}b"

. /software/shared/simumeca/aster/prerequisites/${VERSION}/gcc-mkl-ompi/selena_mpi.sh
