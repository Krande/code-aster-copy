# Code Aster Guideline


## Development
All development dependencies and tasks are handled by `pixi` in a pixi.toml file.

To build, install and test the release version;

```cmd
pixi run build
pixi run install
pixi run test
```

To build and test the debug version you simply add the `d` flag;

```cmd
pixi run buildd
pixi run installd
pixi run testd
```

## MSVC Support

The MSVC support is still experimental.

The MSVC development is to try to leave the original source code as intact as possible.
So as long as this development happens in a separate branch, frequent uses of `git diff` with the git remote `upstream`
should be used to keep track of the changes and ensure minimal disruption to the original source code.

Another goal is to make the code as portable as possible and ultimately make the MSVC version ready for distribution
on conda-forge.

Some changes are necessary to make the code compile with MSVC.
The main additions can be found in the `msvc` directory.

For example, symlinks on Windows require elevated privileges, and is therefore not a viable portable solution
(with conda-forge in mind). Therefore, symlinks are replaced with small c++ redirect modules (`msvc/c_entrypoints`)
which are compiled to dlls and is a replacement for the respective symlink targets on linux and macos.

