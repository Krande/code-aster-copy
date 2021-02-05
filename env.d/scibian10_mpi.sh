# This file set the environment for code_aster.
# Configuration for scibian10 mpi
export WAFBUILD_ENV=$(readlink -n -f ${BASH_SOURCE})

# DEVTOOLS_COMPUTER_ID avoids waf to re-source the environment
export DEVTOOLS_COMPUTER_ID=scibian10

export PREREQ_VERSION=20201218
export PREREQ_PATH=${HOME}/dev/codeaster/prerequisites/${PREREQ_VERSION}/gcc8-ompi2

# force parallel build
export ENABLE_MPI=1

# custom configuration
export CONFIG_PARAMETERS_addmem=2800

export LINKFLAGS="-Wl,--no-as-needed"

# prerequisites paths
export LIBPATH_HDF5="${PREREQ_PATH}/hdf5-1.10.3/lib"
export INCLUDES_HDF5="${PREREQ_PATH}/hdf5-1.10.3/include"
export LD_LIBRARY_PATH="${LIBPATH_HDF5}:${LD_LIBRARY_PATH}"

export LIBPATH_MED="${PREREQ_PATH}/med-4.1.0/lib"
export INCLUDES_MED="${PREREQ_PATH}/med-4.1.0/include"
export PYPATH_MED="${PREREQ_PATH}/med-4.1.0/lib/python3.7/site-packages"
export PATH="${PREREQ_PATH}/med-4.1.0/bin:${PATH}"
export LD_LIBRARY_PATH="${LIBPATH_MED}:${LD_LIBRARY_PATH}"
export PYTHONPATH="${PYPATH_MED}:${PYTHONPATH}"

export LIBPATH_METIS="${PREREQ_PATH}/metis-5.1.0_aster4/lib"
export INCLUDES_METIS="${PREREQ_PATH}/metis-5.1.0_aster4/include"
export LD_LIBRARY_PATH="${LIBPATH_METIS}:${LD_LIBRARY_PATH}"

export TFELHOME="${PREREQ_PATH}/mfront-3.2.1"
export TFELVERS="3.2.1"
export LIBPATH_MFRONT="${PREREQ_PATH}/mfront-3.2.1/lib"
export INCLUDES_MFRONT="${PREREQ_PATH}/mfront-3.2.1/include"
export PYPATH_MFRONT="${PREREQ_PATH}/mfront-3.2.1/lib/python3.7/site-packages"
export PATH="${PREREQ_PATH}/mfront-3.2.1/bin:${PATH}"
export LD_LIBRARY_PATH="${LIBPATH_MFRONT}:${LD_LIBRARY_PATH}"
export PYTHONPATH="${PYPATH_MFRONT}:${PYTHONPATH}"

export PATH="${PREREQ_PATH}/homard-11.12_aster2/bin:${PATH}"

export LIBPATH_SCOTCH="${PREREQ_PATH}/scotch-6.0.4_aster7/lib"
export INCLUDES_SCOTCH="${PREREQ_PATH}/scotch-6.0.4_aster7/include"
export LD_LIBRARY_PATH="${LIBPATH_SCOTCH}:${LD_LIBRARY_PATH}"

export LIBPATH_MUMPS="${PREREQ_PATH}/mumps-5.2.1_consortium_aster3/lib"
export INCLUDES_MUMPS="${PREREQ_PATH}/mumps-5.2.1_consortium_aster3/include"
export LD_LIBRARY_PATH="${LIBPATH_MUMPS}:${LD_LIBRARY_PATH}"

export PATH="${PREREQ_PATH}/miss3d-6.7_aster5/bin:${PATH}"

export LIBPATH_MEDCOUPLING="${PREREQ_PATH}/medcoupling-V9_6_asterxx_0/lib"
export INCLUDES_MEDCOUPLING="${PREREQ_PATH}/medcoupling-V9_6_asterxx_0/include"
export PYPATH_MEDCOUPLING="${PREREQ_PATH}/medcoupling-V9_6_asterxx_0/lib/python3.7/site-packages"
export LD_LIBRARY_PATH="${LIBPATH_MEDCOUPLING}:${LD_LIBRARY_PATH}"
export PYTHONPATH="${PYPATH_MEDCOUPLING}:${PYTHONPATH}"

export PATH="${PREREQ_PATH}/ecrevisse-3.2.2/bin:${PATH}"

export PATH="${PREREQ_PATH}/gmsh-2.12.0-Linux64/bin:${PATH}"

export PATH="${PREREQ_PATH}/grace-0.0.1/bin:${PATH}"

export PYPATH_ASRUN="${PREREQ_PATH}/asrun-2020.0.1/lib/python3.7/site-packages"
export PATH="${PREREQ_PATH}/asrun-2020.0.1/bin:${PATH}"
export PYTHONPATH="${PYPATH_ASRUN}:${PYTHONPATH}"

export PYPATH_MPI4PY="${PREREQ_PATH}/mpi4py-3.0.3/lib/python3.7/site-packages"
export PYTHONPATH="${PYPATH_MPI4PY}:${PYTHONPATH}"

export LIBPATH_PARMETIS="${PREREQ_PATH}/parmetis-4.0.3_aster3/lib"
export INCLUDES_PARMETIS="${PREREQ_PATH}/parmetis-4.0.3_aster3/include"
export LD_LIBRARY_PATH="${LIBPATH_PARMETIS}:${LD_LIBRARY_PATH}"

export LIBPATH_MATH="${PREREQ_PATH}/scalapack-2.1.0/lib"
export LD_LIBRARY_PATH="${LIBPATH_MATH}:${LD_LIBRARY_PATH}"

export LIBPATH_PETSC="${PREREQ_PATH}/petsc-3.12.3_aster/lib"
export INCLUDES_PETSC="${PREREQ_PATH}/petsc-3.12.3_aster/include"
export LD_LIBRARY_PATH="${LIBPATH_PETSC}:${LD_LIBRARY_PATH}"

export LIB_BOOST="boost_python3"
export LIBPATH_BOOST="/usr/lib"
export INCLUDES_BOOST="/usr/include"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${LIBPATH_BOOST}"
