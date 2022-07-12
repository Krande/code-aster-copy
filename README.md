# code_aster source files

[![Documentation Status](https://readthedocs.org/projects/codeaster/badge/?version=latest)](https://codeaster.readthedocs.io/en/latest/?badge=latest)

> Home page: <https://www.code-aster.org/>
>
> Forge: <https://gitlab.com/codeaster/>
>
> code_aster repositories are moving to Git, the work is still in progress...

**code_aster** source files are dispatched into 3 repositories.

- [src][1]: containing Python, C/C++,
  Fortran source files, its build scripts and most of the testcases,
- *validation*: few testcase files with proprietary datas,
- *data*: material datas that can not be freely distributed.

Other independent repositories exist:

- [devtools][2]: contains helper scripts.
- [changelog][3]: publishes the changelog of each incremental version.

## Content of the [src][1] repository

This repository contains to source files of code_aster and its build scripts.

Usually the repositories are cloned under `$HOME/dev/codeaster`.

## Branches and tags

The branches are:

- `main`: The development branch where the current work goes.

- `v15`: The maintenance branch for the stable versions. Only fixes are added,
no new features.

- `v14`: The branch of the old stable version. Not updated anymore.

Each published version is tagged with its number. Examples: 15.4.8, 16.0.9.

Two tags are used aliases and moved when new versions are published:

- `stable`: The latest frozen state of the stable version in the
maintenance branch (ex. 15.5.0).

- `testing`: The latest frozen state of the development version in the
development branch (ex. 16.1.0).

Two names are often used in discussions to identify the code during its
enhancements:

- `unstable`: The head of the development branch.

- `stable-updates`: The head of the maintenance branch.

## Installation

code_aster needs some prerequisites.

Singularity containers are available from the [website][9] in the section Download/salome_meca.
The current version salome_meca and the prerequisites are provided in these containers.
Of course the container must be updated each time that new prerequisites are required
by the development version.

See the [Installation and Development][4] documentation and its changelog page to select
the correct image to be used.

### Example

For the version 16.2.0, the container name is [salome_meca-lgpl-2021.1.0-1-20220405-scibian-9.sif][5].

```bash
$ mkdir $HOME/containers && cd $HOME/containers

# download and the container image (.sif) here
$ singularity run --app install salome_meca-lgpl-2021.1.0-1-20220405-scibian-9.sif

# build code_aster in the container environment
$ cd $HOME/dev/codeaster/src
$ $HOME/containers/salome_meca-lgpl-2021.1.0-1-20220405-scibian-9 --shell

Singularity> ./waf configure install

# or, with the embeeded makefile
Singularity> make bootstrap
```

[1]: ../../../../src
[2]: ../../../../devtools
[3]: ../../../../changelog
[4]: https://gitlab.com/codeaster-opensource-documentation/opensource-installation-development
[5]: https://www.code-aster.org/FICHIERS/singularity/salome_meca-lgpl-2021.1.0-1-20220405-scibian-9.sif
[9]: https://www.code-aster.org/
