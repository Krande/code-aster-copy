#!/bin/bash

_usage()
{
    echo "Check that the documentation is built without error and that automatically"
    echo "generated files are committed."
    echo
    echo "'waf install' (parallel build) must have been run"
    echo "just before. The waf script can be changed but it must point to a"
    echo "parallel configuration."
    echo "This script should be run in the same environment as 'waf install' was."
    echo
    echo "Usage: $(basename $0) [options]"
    echo
    echo "Options:"
    echo
    echo "  --help (-h)            Print this help information and exit."
    echo "  --waf script           Define the script to be used (default: ./waf_mpi)."
    echo "  --use-debug            Use the debug build (default: 'release' using 'waf_mpi')."
    echo "  --verbose (-v)         Show commands output."
    echo
    exit "${1:-1}"
}

check_docs_main()
{
    local waf=./waf_mpi
    local variant="release"
    local verbose=0

    OPTS=$(getopt -o hv --long help,verbose,use-debug,waf: -n $(basename $0) -- "$@")
    if [ $? != 0 ] ; then
        _usage >&2
    fi
    eval set -- "$OPTS"
    while true; do
        case "$1" in
            -h | --help ) _usage ;;
            -v | --verbose ) verbose=1 ;;
            --use-debug ) variant="debug" ;;
            --waf ) waf="$2" ; shift ;;
            -- ) shift; break ;;
            * ) break ;;
        esac
        shift
    done

    if [ -d .git ]; then
        getstatus="git status --porcelain"
    elif [ -d .hg ]; then
        getstatus="hg status -ardm"
    fi
    [ -z "${getstatus}" ] && echo "not a repository" && exit 1

    # with waf_debug, doc and doc_debug are identical
    local suffix=""
    if [ ${variant} = "debug" ]; then
        suffix="_debug"
    fi
    local log=$(mktemp)

    echo -n "Checking docs... "

    (
        (
            printf "\nGenerate html documentation...\n"
            ${waf} doc${suffix} --force-doc
            return $?
        )
        iret=$?
        [ ${iret} -ne 0 ] && return 1

        if [ $(${getstatus} | wc -l) != 0 ]; then
            printf "\nChanges must be committed:\n"
            ${getstatus}
            return 1
        fi

        return 0

    ) >> ${log} 2>&1
    ret=$?

    test "${ret}" = "0" && ok=ok || ok=failed
    echo ${ok}

    if [ "${ok}" != "ok" ] || [ ${verbose} -eq 1 ]
    then
        printf "\nOutput+Error:\n"
        cat ${log}
    fi
    rm -f ${log}

    return ${ret}
}

check_docs_main "$@"
exit $?
