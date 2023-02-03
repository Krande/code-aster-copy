#!/bin/bash -e

usage()
{
    echo "usage: ${0} [arguments]"
    echo
    echo "  This script allows to debug the pipeline locally."
    echo
    echo "arguments:"
    echo "  --all       run all jobs"
    echo "  --prepare   run 'prepare' job"
    echo "  --compile   run 'compile' job"
    echo "  --doc       run 'doc_html' job"
    echo "  --check     run 'check_source' job"
    echo "  --test      run 'minimal_test' job"
    echo
    echo "The repository should be partially cloned (gitlab uses depth=50)."
    echo "From a local repository:"
    echo "  git clone --depth=1 file://path-to-local-clone/src testCI"
    echo "  cd testCI"
    echo "  .gitlabci/debug-ci.sh ..."
}

export CI_PROJECT_DIR=$(pwd)
export CI_SERVER_URL=https://gitlab.pleiade.edf.fr
export CI_REPOSITORY_URL=$(pwd)
export CI_MERGE_REQUEST_SOURCE_BRANCH_NAME=$(git rev-parse --abbrev-ref HEAD)

export DEBUG_CI=1
export ARTF=/tmp
[ -d /local00/tmp ] && export ARTF=/local00/tmp
export HOME=${ARTF}/home
mkdir -p ${HOME}

# variables (: -> =, ' -> ", +export)
export MINIO_URL=https://minio.retd.edf.fr
export SIF=runner.sif
export BUILD=mpi
export GIT_SSL_NO_VERIFY="true"
export DEVTOOLS_ROOT="$CI_PROJECT_DIR/devtools"

SINGULARITY_CMD=(
    "singularity"
    "exec"
    "--cleanenv"
    "--env" "DEVTOOLS_ROOT=${DEVTOOLS_ROOT}"
    "--home=${HOME}"
    "--bind" "$(pwd)"
    "--pwd" "$(pwd)"
    "${SIF}"
)

_cleanup() {
    cd ${CI_PROJECT_DIR}
    echo "+++ cleanup $(pwd)..."
    rm -rf $(git status --porcelain --untracked-files=normal --ignored 2>&1 \
        | egrep '^(\?\?|!!)' | awk '{print $2}')
}

_extract() {
    echo "+++ extracting artifacts from job '$1'..."
    tar xf ${ARTF}/${1}-artifacts.tar
}

_store() {
    echo "+++ creating artifacts for job '$1'..."
    tar cf ${ARTF}/${1}-artifacts.tar $(cat artifacts)
}

do_prepare() {
    _cleanup
    .gitlabci/prepare.sh

    cat << eof > artifacts
${SIF}
devtools/*
eof
    _store prepare
    _cleanup
}

do_compile() {
    _cleanup
    _extract prepare
    "${SINGULARITY_CMD[@]}" .gitlabci/compile.sh

    cat << eof > artifacts
.lock*
build/mpi*/config.log
build/mpi*/c4che/*
build/mpi*/*/*.h
build/mpi*/*/*.py
build/mpi*/*/code_aster/*.py
build/mpi*/*/*/*.so
build/mpi*/*/*.mod
install/*
eof
    _store compile
    _cleanup
}

do_doc_html() {
    _cleanup
    _extract prepare
    _extract compile
    "${SINGULARITY_CMD[@]}" .gitlabci/doc_html.sh

    cat << eof > artifacts
install/share/doc/html/*
eof
    _store doc_html
    _cleanup
}

do_check_source() {
    _cleanup
    _extract prepare
    _extract compile
    "${SINGULARITY_CMD[@]}" .gitlabci/check_source.sh
    _cleanup
}

do_minimal_test() {
    _cleanup
    _extract prepare
    _extract compile
    "${SINGULARITY_CMD[@]}" \
        .gitlabci/test.sh -R "(asrun0|mumps02b|supv002|vocab0|zzzz503n)" -LE need_data --resutest=results_mini

    cat << eof > artifacts
results_mini/run_testcases.xml
results_mini/Testing/Temporary/*
results_mini/*
eof
    _store minimal_test
    _cleanup
}

pipeline() {
    local prepare=0
    local compile=0
    local doc_html=0
    local check_source=0
    local minimal_test=0

    OPTS=$(getopt -o h --long help,all,prepare,compile,doc,check,test,clean -n $(basename $0) -- "$@")
    if [ $? != 0 ] ; then
        _error "invalid arguments." >&2
    fi
    eval set -- "$OPTS"
    while true; do
        case "$1" in
            -h | --help) usage; exit 1 ;;
            --all) prepare=1; compile=1; doc_html=1; check_source=1; minimal_test=1 ;;
            --prepare ) prepare=1 ;;
            --compile ) compile=1 ;;
            --doc ) doc_html=1 ;;
            --check ) check_source=1 ;;
            --test ) minimal_test=1 ;;
            --clean ) _cleanup; exit ;;
            -- ) shift; break ;;
            * ) break ;;
        esac
        shift
    done
    [ $# -ne 0 ] && _error "unexpected argument(s): ${@}"

    if [ "$0" != ".gitlabci/debug-ci.sh" ]; then
        echo "--- must be run in the repository working directory with: .gitlabci/debug-ci.sh"
        exit 1
    fi
    if [ ! -z ${SINGULARITY_NAME} ]; then
        echo "--- $0 must not be run in a singularity container!"
        exit 1
    fi
    if [ $(git status --porcelain -uno | wc -l) != "0" ]; then
        echo "--- there are uncommitted changes!"
        exit 1
    fi

    echo "+++ debugging pipeline: branch=${CI_MERGE_REQUEST_SOURCE_BRANCH_NAME}"
    [ ${prepare} -eq 1 ] && do_prepare
    [ ${compile} -eq 1 ] && do_compile
    [ ${doc_html} -eq 1 ] && do_doc_html
    [ ${check_source} -eq 1 ] && do_check_source
    [ ${minimal_test} -eq 1 ] && do_minimal_test

    echo "+++ run: 'git fetch --unshallow' to complete the repository"
}

pipeline "$@"
exit $?
