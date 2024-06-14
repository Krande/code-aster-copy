# Conda Windows Compilation

You'll need to install MSVC toolchain,

for example by doing `winget install --id=Microsoft.VisualStudio.2022.BuildTools  -e` and add c++ desktop 
development, and specify the MSVC c++ 14.38 SDK for x86 type processors.

In order to compile Code Aster using conda dependencies for windows, you need to perform the following steps:
1. Install conda `call conda\scripts\install_conda.bat`
2. Install Intel Fortran OneAPI 2024.0 compiler `call conda\scripts\install_ifx.bat`
3. Open the miniforge terminal, cd to this directory and run:

```cmd
mamba env update -f env.debug.yml
```

3. Install VS Build Tools (or the full installation of) VS2022 (https://aka.ms/vs/17/release/vs_BuildTools.exe or https://visualstudio.microsoft.com/downloads/)
3. Create an `.env` file containing the paths to the installation of python, VS Build Tools and Intel Fortran OneAPI
   2024.0

Example `.env` file:

```
CONDA_ROOT=C:\work\miniforge3
INTEL_VARS_PATH=C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env
# BUILD TOOLS PATH
VS_VARS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build
# VS2022 Professional PATH
#VS_VARS_PATH=C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build
PYTHON_ENV_NAME=ca-debug
```

4. Run the batch file `conda_build.bat` to compile Code Aster for Windows

```cmd
conda_build.bat
```

## Windows Dependencies Compiler versions (conda-forge)

Conda forge currently does not support Intel OneAPI fortran. Consequently, focus is shifted to adding support
for LLVM Flang (which is supported on conda-forge). 

| Dependency  | version  | C compiler | C++ compiler | Fortran compiler                 | 
|-------------|----------|------------|--------------|----------------------------------|
| HDF5        | 1.14.4.2 | VS2022     | VS2022       | Intel OneAPI Fortran 2024.0 (^1) |
| MED         | 4.1.0    | VS2022     | VS2022       | Intel OneAPI Fortran 2024.0 (^1) |
| MEDCOUPLING | 9.10.0   | VS2022     | VS2022       | N/A                              |
| MFront      | 4.2.0    | CLANG-CL   | CLANG-CL     | LLVM Flang                       |
| MGIS        | 2.2.0    | VS2019     | VS2019       | LLVM Flang                       |
| METIS       | 5.1.0    | VS2022     | VS2022       | N/A                              |
| SCOTCH      | 7.0.4    | VS2022     | VS2022       | LLVM Flang                       |
| MUMPS       | 5.7.0    | VS2022     | VS2022       | LLVM Flang                       |
| Code Aster  | 17.0.10  | VS2022     | VS2022       | Intel OneAPI Fortran 2024.0 (^1) |

^1: Awaiting LLVM Flang fix -> https://github.com/llvm/llvm-project/issues/89403 

## OpenMP error related to calculation of `omp_get_max_threads`

If you encounter a large spike in memory usage, it is likely caused by the calculation of `omp_get_max_threads` in
`bibfor/supervis/superv_module.F90`. If you look at `Nombre de processus OpenMP utilisés : 1` in the Code Aster
output, you can see that the number of threads is set to 1. However, when a large spike of memory is observed,
the number of OpenMP threads is likely a much higher number. 

The actual culprit for this error is still not 100% clear, but it is likely related to the calculation of the number
of threads in the `bibfor/supervis/superv_module.F90` file.

Originally I thought I had solved this issue for good when I changed MUMPS to use MKL64 and enable OpenMP, but that
might not be the case.

It might be related to the PATH environment variable and the Intel oneAPI vars set globally. Removing the paths and 
restarting the computer seemed to have worked last time I encountered this issue.




```
                       -- CODE_ASTER -- VERSION : DÉVELOPPEMENT (unstable) --                       
                               Version 17.0.99 modifiée le 07/06/2024                               
                                    révision n/a - branche 'n/a'                                    
                                   Copyright EDF R&D 1991 - 2024                                    
                                                                                                    
                              Exécution du : Sat Jun  8 10:02:55 2024                               
                                        Nom de la machine :                                         
                                          DESKTOP-NJ4G2LS                                           
                                        Architecture : 64bit                                        
                                     Type de processeur : AMD64                                     
                                      Système d'exploitation :                                      
                                     Windows-10-10.0.22631-SP0                                      
                                 Langue des messages : nb (cp1252)                                  
                                     Version de Python : 3.11.9                                     
                                     Version de NumPy : 1.23.5                                      
                                     Parallélisme MPI : inactif                                     
                                    Parallélisme OpenMP : actif                                     
                              Nombre de processus OpenMP utilisés : 1                               
                               Version de la librairie HDF5 : 1.14.4                                
                                Version de la librairie MED : 4.1.1                                 
                               Version de la librairie MFront : 4.2.0                               
                               Version de la librairie MUMPS : 5.6.2                                
                                  Librairie PETSc : non disponible                                  
                               Version de la librairie SCOTCH : 7.0.4  
```

## Debugging

The conda directory contains a series of manual lib and link scripts for testing and debugging the creation of 
code aster DLL's. 

Use the `--bibc` flag to re-link the compiled .o sources for bibc.

## Resources

https://github.com/scipy/scipy/pull/16957