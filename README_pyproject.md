# Dependencies (Experimental)

**!!This feature is still experimental!!**

Dependencies are defined in pyproject.toml file

There are two kinds of dependencies :

- pypi which are defined in standard pyproject section according to the PEPs.
Those dependencies can be installed either with [uv astral](https://docs.astral.sh/uv/)
or [pixi](https://pixi.sh/latest/)
- conda (python or other) which can only be installed with [pixi](https://pixi.sh/latest/)

# uv vs pixi

uv only handles pypi (pure python) dependencies.
pixi uses uv for pypi dependencies and conda for others.
pixi implements tasks using Deno shell (lightweight posix shell compatible with linux windows and macos).

- if you only want to use python dependencies uv is just fine (for example for autocompletion in your editor)
- if you want requirements installations for aster compilation + tasks you should pixi

# Pixi Tasks (Experimental)

**!!This feature is still experimental!!**

We provide several tasks to ease compilation, when executed those tasks will update environments
and then run tasks.
Tasks can depend on other ones for instance test task depends on build -> build will be triggered
every time you run test.

[pixi tasks](https://pixi.sh/latest/workspace/advanced_tasks/)

- To list all available tasks :

```bash
pixi task list
```

> [!NOTE]
> `pixi r` is a shortcut for `pixi run`

As a starting point, here is how you can install aster using pixi tasks:

```bash
pixi r configure
pixi r install
```

Shortcut to run_aster executable:

```bash
pixi r start ....
```

or if you want to install aster and then use run_ctest (usefull for development iterations):

```bash
pixi r test -R "zzzz" ....
```

# Pixi Environments

Pixi provides the capacity to handle several environments in a single project.
It's usefull to isolate dependencies and avoid dependency conflicts.

To list all available project pixi environments:

```bash
pixi workspace environment list
```

[pixi environments](https://pixi.sh/latest/workspace/environment/)

# EDF specific

Edf users can install pixi using this procedure [pixi install edf](https://gitlab.pleiade.edf.fr/codeaster/lab/pixi-install-helper/-/blob/main/README.md?ref_type=heads)

