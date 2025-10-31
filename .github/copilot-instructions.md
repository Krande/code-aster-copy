# Code Aster Guideline

You are an expert in mixed C/C++/Fortran/Python code, the MSVC toolchain and the WAF build system.
The agent shell is powershell when using windows.

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

The `build` and `buildd` on windows call the `msvc/scripts/conda_manual_build.bat`  
which can take different arguments.

* `--use-log`: Creates a log file for the installation phase (the configure phase is still logged to the console)

All build commands are stored in `build\int64\debug\compile_commands.json` (can be debug/release depending on the build)

## MSVC Support

### Code Divergence
The MSVC development is to trying to leave the original source code as intact as possible.
So as long as the development happens in a separate branch, frequent uses of `git diff` with the git remote `upstream`
should be used to continuously monitor the divergence from the main branch and to always ensure minimal disruption
to the original source code.

As an attempt to isolate changes to the main source code MSVC configuration and additional files
are placed in the `msvc` directory.

### Portability
Another goal is to make the code as portable as possible and ultimately make the MSVC version ready for distribution
on conda-forge.

So one of the main issues of portability in porting Code_Aster to MSVC and prepare it for conda-forge
is the use of symlinks between the modules when compiled using shared libraries.

Symlinks on Windows require elevated privileges, and is therefore not a viable portable solution
(with conda-forge in mind). Therefore, symlinks are replaced with small c++ redirect modules (`msvc/c_entrypoints`)
which are compiled to the respective symlink targets used on linux or MinGW/MSYS2.

### Symbol exports

Symbol exports. The MSVC toolchain is a bit different from the Linux toolchain, and therefore
the symbol exports are not as straightforward as on Linux.

In `msvc` directory there are scripts for producing the `*.def` files necesary to export symbols on Windows.

Pixi tasks are available to recreate the def files if needed:

```cmd
pixi run def-bibfor
```


```cmd
pixi run def-bibcpp
```