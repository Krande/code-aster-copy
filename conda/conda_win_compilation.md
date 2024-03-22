# Conda Windows Compilation

In order to compile Code Aster using conda dependencies for windows, you need to perform the following steps:

1. Open a terminal, cd to this directory and run:

```cmd
mamba env update -f environment.yml
```

2. Install VS Build Tools (or the full installation of) VS2022 and Intel Fortran OneAPI 2024.0
3. Create an `.env` file containing the paths to the installation of python, VS Build Tools and Intel Fortran OneAPI 2024.0

Example `.env` file:
```
CONDA_ROOT=C:\work\miniforge3
INTEL_VARS_PATH=C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env
VS_VARS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build
```

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


## Packages


### Intel Fortran
```
# wget and 
curl https://registrationcenter-download.intel.com/akdlm/IRC_NAS/3a64aab4-3c35-40ba-bc9c-f80f136a8005/w_fortran-compiler_p_2024.0.2.27_offline.exe -o w_fortran-compiler_p_2024.0.2.27_offline.exe
w_fortran-compiler_p_2024.0.2.27_offline.exe -s

```