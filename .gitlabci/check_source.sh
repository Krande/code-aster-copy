#!/bin/bash -e

if [ ! -z "${GITLAB_CI}" ]; then
    echo "+ fetching '${REFREV}' branch..."
    git branch -D ${REFREV} || true
    git fetch --depth=50 origin ${REFREV}
    git branch ${REFREV} FETCH_HEAD
fi
base=$(git merge-base ${REFREV} HEAD)

echo "+ printing all branches..."
git branch -av

echo "+ checking changes..."
git diff --name-only ${base}

echo "+ calling aslint..."
./devtools/bin/aslint --force --reponame=${CI_PROJECT_NAME} --repo ${base}

flist=list_issues.txt
fmsg=$(mktemp tmp.msg.XXXXXXXX)
trap "rm -f ${fmsg}" EXIT

echo "+ checking commit messages..."
for rev in $(git rev-list ${base}..HEAD)
do
    echo "- checking revision ${rev}..."
    git log --format=%s -1 ${rev} > ${fmsg}
    grep -Eo '\[.*\]' ${fmsg} | grep -Eo '#[0-9]+' | sed -e 's/#//' >> ${flist}
    ./devtools/bin/hooks/git_hook --commit-msg ${fmsg}
done
