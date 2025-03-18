# This file set the environment for code_aster.
# Configuration to use prerequisites from salome_meca installation

if [ ! -z "${ROOT_SALOME_INSTALL}" ]; then
    export DEVTOOLS_COMPUTER_ID=salome_meca
    export ASTER_BUILD=std
    export ENABLE_MPI=0

    libpaths=$( sed -e 's/:/ /g' <<< "${LD_LIBRARY_PATH}" )
    components=( "HDF5" "MED:MEDFILE" "MFRONT:MFRONT_TESTING" "MGIS:MGIS_TESTING"
                 "METIS:METIS_ASTER" "SCOTCH:SCOTCH_TESTING"
                 "MUMPS:MUMPS_TESTING" )
    # "PARMETIS" "PETSC" "PETSC4PY" "NUMPY" "MPI"
    for comp in ${components[@]}
    do
        prod=$(awk -F: '{if (NF==2) {print $2} else {print $1}}' <<< "${comp}")
        comp=$(awk -F: '{print $1}' <<< "${comp}")
        rootdir="$(eval echo \$${prod}_ROOT_DIR)"
        if [ "${comp}" = "MUMPS" ] || [ "${comp}" = "SCOTCH" ]; then
            rootdir="${rootdir}/SEQ"
        fi
        export LIBPATH_${comp}="${rootdir}/lib"
        export INCLUDES_${comp}="${rootdir}/include"
        if [ -z "${INCLUDES}" ]; then
            export INCLUDES="${rootdir}/include"
        else
            export INCLUDES="${INCLUDES}:${rootdir}/include"
        fi
        if [ "${comp}" = "MUMPS" ]; then
            export INCLUDES_MUMPS="${INCLUDES_MUMPS} ${rootdir}/include_seq"
            export INCLUDES="${INCLUDES}:${rootdir}/include_seq"
        fi
    done
    export TFELHOME="${MFRONT_TESTING_ROOT_DIR}"
    export TFELVERS="4.2.2"
else
    echo "SALOME environment not found, please use 'salome shell' before starting configure."
fi
