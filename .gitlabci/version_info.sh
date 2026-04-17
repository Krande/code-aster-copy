#!/bin/bash
set -e

TAG=$(git describe --tags --abbrev=0 2>/dev/null || echo "no-tag")
REVISION=$(git rev-parse HEAD)
DISTANCE=$(git rev-list "${TAG}"..HEAD --count 2>/dev/null || echo "0")
DATE=$(git show -s --format=%cd --date=format:%d/%m/%Y HEAD)
BRANCH=$(git branch --show-current || echo "detached")

echo "Version info:"
echo "($TAG, $REVISION, $BRANCH, $DATE, $DISTANCE)"

# GitLab dotenv
cat <<EOF > version.env
VERSION_TAG=$TAG
VERSION_REVISION=$REVISION
VERSION_BRANCH=$BRANCH
VERSION_DATE=$DATE
VERSION_DISTANCE=$DISTANCE
EOF
