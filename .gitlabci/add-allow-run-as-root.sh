#!/bin/bash

sed -i -e "s/mpiexec *: *mpiexec/mpiexec: mpiexec --allow-run-as-root/g" \
    install/share/aster/config.yaml
