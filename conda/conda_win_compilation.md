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

## Windows Dependencies Compiler versions

| Dependency  | version | C compiler | C++ compiler | Fortran compiler    | 
|-------------|---------|------------|--------------|---------------------|
| HDF5        | 1.10.7  | VS2022     | VS2022       | N/A                 |
| MED         | 4.1.0   | VS2022     | VS2022       | N/A                 |
| MEDCOUPLING | 9.10.0  | VS2022     | VS2022       | N/A                 |
| MFront      | 4.2.0   | CLANG-CL   | CLANG-CL     | N/A                 |
| MGIS        | 2.2.0   | VS2022     | VS2022       | N/A                 |
| METIS       | 5.1.0   | VS2022     | VS2022       | N/A                 |
| SCOTCH      | 7.0.4   | VS2022     | VS2022       | N/A                 |
| Code Aster  | 17.0.10 | VS2022     | VS2022       | Intel OneAPI 2024.0 |


