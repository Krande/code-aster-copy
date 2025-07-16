#!/bin/bash -e

flist=list_issues.txt

echo "+ checking that expected documents have been committed..."
./devtools/bin/maint/check_expected_documents ${flist}

echo "+ checking issues status..."
./devtools/bin/maint/check_issue_status --expected=valide_EDA ${flist}
