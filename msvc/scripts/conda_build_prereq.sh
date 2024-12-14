#!/bin/bash
set -e

# Set the name of the current prerequisites tarball
#fname=codeaster-prerequisites-20240327-oss
fname=codeaster-prerequisites-20221225-oss

fp=${fname}.tar.gz
root_url=https://github.com/Krande/condapackaging/releases/download
url_tar=$root_url/$fname/$fname.tar.gz

# Make a temporary subdirectory "temp/packages" which is also set as a env variable
mkdir -p temp/packages
export PACKAGES_DIR=$(pwd)/temp/packages
export PREREQUISITES_DIR=$PACKAGES_DIR/$fname

# If file does not exists, download it to the PACKAGES_DIR
if [ ! -f $PACKAGES_DIR/$fp ]; then
  wget -O $PACKAGES_DIR/$fp $url_tar
fi

# Extract the tarball if not already extracted
if [ ! -d $PREREQUISITES_DIR ]; then
  tar -xzf $PACKAGES_DIR/$fp -C $PACKAGES_DIR
fi

