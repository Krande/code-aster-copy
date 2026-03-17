# Build with `spin` commands (Experimental)

**!!This feature is still experimental!!**

`spin` commands are available to build code_aster.

> _The goal of spin is to provide a simple, user-friendly, extendable interface for common development tasks._
>
> More on [Scientific Python INcantations (spin)](https://github.com/scientific-python/spin).

**NB:** `spin` is not provided within code_aster. It must be installed (in a venv using `pip` for example).


### Examples

- Basic example:

```shell
$ spin configure
$ spin install
$ spin doc
$ spin test ssnv128p
```

- Options are passed to underlying commands:

```shell
$ spin configure --prefix=/opt/aster
$ spin install -j 16 --install-tests --fast
$ spin test mumps01a mumps02a
```

- All in one command (bootstrap):

```shell
$ spin bootstrap
```

- Available commands:

```shell
$ spin
Usage: spin [OPTIONS] COMMAND [ARGS]...

  Developer tool for codeaster-src

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  configure  Configure the project
  build      Build the project
  install    Build and install the project
  doc        Build the embedded documentation
  test       Run one or more testcases
  bootstrap  Run all steps: configure, install and doc
  distclean  Run distclean.
```
