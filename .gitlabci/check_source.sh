#!/bin/bash -e

echo "+ fetching 'main' branch..."
git branch -D main || true
git fetch --depth=1 origin main
git branch main FETCH_HEAD

echo "+ printing all branches..."
git branch -av

echo "+ checking changes..."
git diff --name-only main

echo "+ calling aslint..."
./devtools/bin/aslint --force --reponame=src --repo main
