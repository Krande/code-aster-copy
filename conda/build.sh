#!/bin/bash

export CLICOLOR_FORCE=1

echo "mpi = ${mpi}"
echo "mpi_prefix = ${mpi_prefix}"
echo "build_type = ${build_type}"
mpi=${mpi_prefix}

export CONFIG_PARAMETERS_addmem=2000
export TFELHOME=${PREFIX}

export LIBPATH_METIS="${PREFIX}/lib"
export INCLUDES_METIS="${PREFIX}/include"

export LIBPATH_PETSC="${PREFIX}/lib"
export INCLUDES_PETSC="${PREFIX}/include"

export INCLUDES_BOOST=${PREFIX}/include
export LIBPATH_BOOST=${PREFIX}/lib
export LIB_BOOST="libboost_python$CONDA_PY"

export INCLUDES_MUMPS="${PREFIX}/include"
export LIBPATH_MUMPS="${PREFIX}/lib"

export INCLUDES_MED="${PREFIX}/include"
export LIBPATH_MED="${PREFIX}/lib"

export LIBPATH_MEDCOUPLING="${PREFIX}/lib"
export INCLUDES_MEDCOUPLING="${PREFIX}/include"
export PYPATH_MEDCOUPLING=${SP_DIR}

python conda/scripts/update_version.py

mpi_type=std
if [[ "$mpi" != "nompi" ]]; then
  mpi_type=mpi
fi

if [[ "${build_type}" == "debug" ]]; then
    echo "Debugging Enabled"
    export CFLAGS="-g -O0 ${CFLAGS}"
    export CXXFLAGS="-g -O0 ${CXXFLAGS}"
    export FCFLAGS="-g -O0 ${FCFLAGS}"
    build_type=debug
else
    echo "Debugging Disabled"
fi

export FCFLAGS="-fdefault-integer-8 ${FCFLAGS}"
export FFLAGS="-fdefault-integer-8 ${FFLAGS}"

# if gfortran version > 9, we need to conditionally add -fallow-argument-mismatch
# to avoid mismatch errors related to floats and integer types
major_version=$($FC -dumpversion | awk -F. '{print $1}')
if [[ $major_version -gt 9 ]]; then
  echo "adding -fallow-argument-mismatch to FCFLAGS"

  export FCFLAGS="-fallow-argument-mismatch ${FCFLAGS}"
else
  echo "FCFLAGS: $FCFLAGS"
fi


if [[ "$mpi" == "nompi" ]]; then
  # Install for standard sequential
  waf \
    --use-config-dir=${SRC_DIR}/config/ \
    --prefix="${PREFIX}" \
    --med-libs="med medC medfwrap medimport" \
    --enable-med \
    --enable-hdf5 \
    --enable-mumps \
    --enable-metis \
    --enable-scotch \
    --enable-mfront \
    --libdir="${PREFIX}/lib" \
    --install-tests \
    --disable-mpi \
    --disable-petsc \
    --without-hg \
    configure

  if [[ "${build_type}" == "debug" ]]; then
      waf install_debug -v
  else
      echo "Debugging Disabled"
      waf install
  fi
else
  export PYTHONPATH="$PYTHONPATH:${PREFIX}/lib"
  export CONFIG_PARAMETERS_addmem=4096

  export ENABLE_MPI=1
  export CC=mpicc
  export CXX=mpicxx
  export FC=mpif90
  export F77=mpif77
  export F90=mpif90
  export OPAL_PREFIX=${PREFIX}

  waf configure \
    --use-config-dir=${SRC_DIR}/config/ \
    --enable-med \
    --enable-hdf5 \
    --enable-mumps \
    --enable-metis \
    --enable-scotch \
    --enable-mfront \
    --med-libs="med medC medfwrap medimport" \
    --prefix="${PREFIX}" \
    --enable-mpi \
    --libdir="${PREFIX}/lib" \
    --install-tests \
    --without-hg


  if [[ "${build_type}" == "debug" ]]; then
      waf install_debug
  else
      waf install
  fi
fi

echo "Compilation complete"

export LD_LIBRARY_PATH="${PREFIX}/lib/aster"

# This is for reducing reliance on conda activation scripts.
mv "${PREFIX}/lib/aster/code_aster" "${SP_DIR}/code_aster"
mv "${PREFIX}/lib/aster/run_aster" "${SP_DIR}/run_aster"

# move aster_pkginfo.py and aster_version.py
mv ${SRC_DIR}/build/${build_type}/code_aster/*.py "${SP_DIR}/code_aster/Utilities/"
# note to self. aster.so is symlinked to libaster.so
mv ${PREFIX}/lib/aster/libb*.so "${PREFIX}/lib/"
mv ${PREFIX}/lib/aster/libAsterMFrOff*.so "${PREFIX}/lib/"

mv "${PREFIX}/lib/aster/med_aster.so" "${SP_DIR}/"
mv ${PREFIX}/lib/aster/*.so "${SP_DIR}/"
