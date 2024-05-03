# Conda Windows Compilation

In order to compile Code Aster using conda dependencies for windows, you need to perform the following steps:
1. Install conda `call conda\install_conda.bat`
2. Install Intel Fortran OneAPI 2024.0 compiler `call conda\install_ifx.bat`
3. Open the miniforge terminal, cd to this directory and run:

```cmd
mamba env update -f environment.yml
```

3. Install VS Build Tools (or the full installation of) VS2022 (https://aka.ms/vs/17/release/vs_BuildTools.exe or https://visualstudio.microsoft.com/downloads/)
3. Create an `.env` file containing the paths to the installation of python, VS Build Tools and Intel Fortran OneAPI
   2024.0

Example `.env` file:

```
CONDA_ROOT=C:\work\miniforge3
INTEL_VARS_PATH=C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env
VS_VARS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build
PYTHON_ENV_NAME=codeaster-deps
```

4. Run the batch file `conda_build.bat` to compile Code Aster for Windows

```cmd
conda_build.bat
```

## Windows Dependencies Compiler versions (conda-forge)

Conda forge currently does not support Intel OneAPI fortran. Consequently, focus is shifted to adding support
for LLVM Flang (which is supported on conda-forge). 

| Dependency  | version | C compiler | C++ compiler | Fortran compiler                 | 
|-------------|---------|------------|--------------|----------------------------------|
| HDF5        | 1.10.7  | VS2022     | VS2022       | Intel OneAPI Fortran 2024.0 (^1) |
| MED         | 4.1.0   | VS2022     | VS2022       | Intel OneAPI Fortran 2024.0 (^1) |
| MEDCOUPLING | 9.10.0  | VS2022     | VS2022       | N/A                              |
| MFront      | 4.2.0   | CLANG-CL   | CLANG-CL     | LLVM Flang                       |
| MGIS        | 2.2.0   | VS2022     | VS2022       | LLVM Flang                       |
| METIS       | 5.1.0   | VS2022     | VS2022       | N/A                              |
| SCOTCH      | 7.0.4   | VS2022     | VS2022       | LLVM Flang                       |
| MUMPS       | 5.7.0   | VS2022     | VS2022       | LLVM Flang                       |
| Code Aster  | 17.0.10 | VS2022     | VS2022       | Intel OneAPI Fortran 2024.0 (^1) |

^1: Awaiting LLVM Flang fix -> https://github.com/llvm/llvm-project/issues/89403 


## Debugging

The conda directory contains a series of manual lib and link scripts for testing and debugging the creation of 
code aster DLL's. 

Use the `--bibc` flag to re-link the compiled .o sources for bibc.

## Resources

https://github.com/scipy/scipy/pull/16957