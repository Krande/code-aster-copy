#!/bin/bash -e

refrev=main

echo "+ fetching '${refrev}' branch..."
git branch -D ${refrev} || true
git fetch --depth=1 origin ${refrev}
git branch ${refrev} FETCH_HEAD

echo "+ printing all branches..."
git branch -av

echo "+ checking changes..."
git diff --name-only ${refrev}

echo "+ calling aslint..."
./devtools/bin/aslint --force --reponame=${CI_PROJECT_NAME} --repo ${refrev}

flist=$(mktemp tmp.list.XXXXXXXX)
fmsg=$(mktemp tmp.msg.XXXXXXXX)
trap "rm -f ${fmsg} ${flist}" EXIT

echo "+ checking commit messages..."
for rev in $(git rev-list ${refrev}..HEAD)
do
    echo "- checking revision ${rev}..."
    git log --format=%s -1 ${rev} > ${fmsg}
    grep -Eo '\[.*\]' ${fmsg} | grep -Eo '#[0-9]+' | sed -e 's/#//' >> ${flist}
    ./devtools/bin/hooks/git_hook --commit-msg ${fmsg}
done

echo "+ checking issues status..."
./devtools/bin/maint/check_issue_status --expected=valide_EDA ${flist}

echo "+ checking that expected documents have been committed..."
./devtools/bin/maint/check_expected_documents ${flist}
