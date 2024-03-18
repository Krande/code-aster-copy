# Conda Windows Compilation

In order to compile Code Aster using conda dependencies for windows, you need to perform the following steps:

1. Open a terminal, cd to this directory and run:

```cmd
mamba env update -f environment.yml
```

2. Install VS Build Tools (or the full installation of) VS2022 and Intel Fortran OneAPI 2024.0
3. Update the paths relevant to the installation of the VS Build Tools and Intel Fortran OneAPI 2024.0 in the `conda_build.bat` file
4. Run the batch file `conda_build.bat` to compile Code Aster for Windows

```cmd
conda_build.bat
```
