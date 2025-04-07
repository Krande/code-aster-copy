#!/bin/bash -e

usage()
{
    cmd="${0}"
    [ ! -z "${submit}" ] && cmd="check_ci.sh"
    echo "usage: ${cmd} [arguments]"
    echo
    echo "  This script allows to debug the pipeline locally."
    echo
    if [ ! -z "${submit}" ]; then
        echo "  By default, only '--check', '--doc', '--test' are enabled"
        echo "  because '--compile' and '--mini' change your build directory."
        echo
    fi
    echo "arguments:"
    if [ -z "${submit}" ]; then
        echo "  --all       run all jobs"
        echo "  --prepare   run 'prepare' job"
    fi
    echo "  --compile   run 'compile' job"
    echo "  --mini      run 'minimal_build' job"
    echo "  --check     run 'check_source' and 'check_issues' jobs"
    echo "  --doc       run 'doc_html' job"
    echo "  --test      run 'minimal_test' job"
    echo
    if [ -z "${submit}" ]; then
        echo "The repository should be partially cloned (gitlab uses depth=50)."
        echo "From a local repository:"
        echo "  git clone --depth=1 file://path-to-local-clone/src testCI"
        echo "  cd testCI"
        echo "  .gitlabci/debug-ci.sh ..."
    fi
}

submit="${SUBMIT_MODE}"

export CI_PROJECT_DIR=$(pwd)
export CI_SERVER_URL=https://gitlab.pleiade.edf.fr
export CI_PROJECT_URL=$(pwd)
export CI_COMMIT_REF_NAME=$(git rev-parse --abbrev-ref HEAD)
export REFREV=main

export DEBUG_CI=1
export ARTF=/tmp
[ -d /local00/tmp ] && export ARTF=/local00/tmp
if [ -z "${submit}" ]; then
    export ORIG_HOME=${HOME}
    export HOME=${ARTF}/home
    mkdir -p ${HOME}
fi

# variables (: -> =, ' -> ", +export)
export MINIO_URL=https://minio.retd.edf.fr
export SIF=runner.sif
export ASTER_BUILD=mpi
export GIT_SSL_NO_VERIFY="true"

SINGULARITY_CMD=(
    "singularity"
    "exec"
    "--cleanenv"
    "--home=${HOME}"
    "--bind" "$(pwd)"
    "--pwd" "$(pwd)"
    "${SIF}"
)

_cleanup() {
    [ ! -z "${submit}" ] && return
    cd ${CI_PROJECT_DIR}
    echo "+++ cleanup $(pwd)..."
    rm -rf $(git status --porcelain --untracked-files=normal --ignored 2>&1 \
        | egrep '^(\?\?|!!)' | awk '{print $2}')
}

_extract() {
    [ ! -z "${submit}" ] && return
    echo "+++ extracting artifacts from job '$1'..."
    tar xf ${ARTF}/${1}-artifacts.tar
}

_store() {
    [ ! -z "${submit}" ] && return
    echo "+++ creating artifacts for job '$1'..."
    tar cf ${ARTF}/${1}-artifacts.tar $(cat artifacts)
}

do_prepare() {
    [ ! -z "${submit}" ] && return
    _cleanup
    .gitlabci/prepare.sh

    cat << eof > artifacts
${SIF}
devtools/*
data-src/*
eof
    _store prepare
    _cleanup
}

do_compile() {
    [ ! -z "${submit}" ] && return
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
build/mpi*/*/catalo/elem.1
build/mpi*/*/*.mod
install/*
eof
    _store compile
    _cleanup
}

do_compile_mini() {
    _cleanup
    _extract prepare
    "${SINGULARITY_CMD[@]}" .gitlabci/compile_minimal.sh
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
    "${SINGULARITY_CMD[@]}" .gitlabci/check_issues.sh
    _cleanup
}

do_minimal_test() {
    _cleanup
    _extract prepare
    _extract compile
    "${SINGULARITY_CMD[@]}" \
        .gitlabci/test.sh -R "(asrun0|mumps02b|supv002|vocab0|zzzz503n)" --resutest=results_mini

    cat << eof > artifacts
results_mini/run_testcases.xml
results_mini/Testing/Temporary/*
results_mini/*
eof
    _store minimal_test
    _cleanup
}

check_devtools() {
    mark="${HOME}/.config/aster/devtools_fetch.mark"
    ftest=$(mktemp tmp.devtools_fetch.XXXXXXXX)
    touch -d "-10 days" "${ftest}"
    if [ ! -f "${mark}" ] || [ "${mark}" -ot "${ftest}" ]; then
        echo "+ devtools should be regularly updated"
        printf "do you want to updated 'devtools' now (y/[n])? "
        read answ
        if [ "${answ}" = "y" ]; then
            date > "${mark}"
            (cd devtools && git checkout main && git pull)
        fi
    fi
    rm -f "${ftest}"
}

check_install_env() {
    mark="${HOME}/.config/aster/install_env.mark"
    ftest=$(mktemp tmp.install_env.XXXXXXXX)
    touch -d "-20 days" "${ftest}"
    if [ ! -f "${mark}" ] || [ "${mark}" -ot "${ftest}" ]; then
        echo "+ install_env should be regularly run"
        printf "do you want to run 'install_env' now (y/[n])? "
        read answ
        if [ "${answ}" = "y" ]; then
            date > "${mark}"
            ./devtools/bin/install_env --no-build
        fi
    fi
    rm -f "${ftest}"
}

pipeline() {
    local prepare=0
    local compile=0
    local mini=0
    local doc_html=0
    local check_source=0
    local minimal_test=0
    local native=0
    local inall=1
    [ ! -z "${submit}" ] && inall=0

    OPTS=$(getopt -o h --long help,all,prepare,compile,mini,doc,check,test,clean,native -n $(basename $0) -- "$@")
    if [ $? != 0 ] ; then
        _error "invalid arguments." >&2
    fi
    eval set -- "$OPTS"
    while true; do
        case "$1" in
            -h | --help) usage; exit 1 ;;
            --all) prepare=${inall}; compile=${inall}; mini=${inall}; doc_html=1; check_source=1; minimal_test=1 ;;
            --prepare ) prepare=1 ;;
            --compile ) compile=1 ;;
            --mini ) mini=1 ;;
            --doc ) doc_html=1 ;;
            --check ) check_source=1 ;;
            --test ) minimal_test=1 ;;
            --clean ) _cleanup; exit ;;
            --native ) native=1 ;;
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
        echo "+ running inside a singularity container, 'native' mode enabled"
        native=1
    fi
    if [ -f /.dockerenv ]; then
        echo "+ running inside a docker container, 'native' mode enabled"
        native=1
    fi
    if [ ! -z "${submit}" ]; then
        echo "+ running in submit mode"
        cleanup="${cleanup} echo cleanup... ; rm -rf artifacts list_issues.txt results_mini ;"
        trap "${cleanup}" EXIT
        export GIT_PAGER=cat
        check_devtools
        check_install_env
    fi
    if [ $(git status --porcelain -uno | wc -l) != "0" ]; then
        echo "--- there are uncommitted changes, you should commit or stash them!"
        printf "do you want to continue (y/[n])? "
        read answ
        [ "${answ}" != "y" ] && exit 1
    fi
    if [ ${native} -eq 1 ]; then
        export SIF=""
        SINGULARITY_CMD=()
    fi

    echo "++ running on branch=${CI_COMMIT_REF_NAME}"
    [ ${prepare} -eq 1 ] && do_prepare
    [ ${compile} -eq 1 ] && do_compile
    [ ${mini} -eq 1 ] && do_compile_mini
    [ ${check_source} -eq 1 ] && do_check_source
    [ ${doc_html} -eq 1 ] && do_doc_html
    [ ${minimal_test} -eq 1 ] && do_minimal_test

    if [ -z "${submit}" ]; then
        echo "++ run: 'git fetch --unshallow' to complete the repository"
    fi
}

pipeline "$@"
exit $?
