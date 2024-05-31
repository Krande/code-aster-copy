# Issue Tracker

## 24.05.31: Issue with Windows Compatibility

When compiling the code on Windows, the following error message is displayed:

```
11:53:28 runner 'C:\\work\\Miniforge3\\envs\\codeaster-deps\\python.exe C:\\work\\code\\code-aster-src\\build\\std\\debug\\catalo\\fort.1 --memory=5120'
Traceback (most recent call last):
  File "C:\work\code\code-aster-src\build\std\debug\catalo\fort.1", line 1, in <module>
    from code_aster.Commands import DEBUT, MAJ_CATA, FIN
  File "C:\work\miniforge3\envs\codeaster-deps\Library\lib\aster\code_aster\Commands.py", line 29, in <module>
    from .Utilities.rc import rc
  File "C:\work\miniforge3\envs\codeaster-deps\Library\lib\aster\code_aster\Utilities\__init__.py", line 59, in <module>
    from .ExecutionParameter import ExecutionParameter
  File "C:\work\miniforge3\envs\codeaster-deps\Library\lib\aster\code_aster\Utilities\ExecutionParameter.py", line 50, in <module>
    import libaster
ModuleNotFoundError: No module named 'libaster'
```

This is because symlinking is currently not working. 

By using the following command (using a terminal with admin privileges), the issue can be resolved:

```
python conda\symlink.py --symlink
```

When I re-run the compilation (make sure to disable `waf distclean` in the conda_manual_build.bat), 
the error message is then this.

```
11:59:27 runner 'C:\\work\\Miniforge3\\envs\\codeaster-deps\\python.exe C:\\work\\code\\code-aster-src\\build\\std\\debug\\catalo\\fort.1 --memory=5120'
Traceback (most recent call last):
  File "C:\work\code\code-aster-src\build\std\debug\catalo\fort.1", line 1, in <module>
    from code_aster.Commands import DEBUT, MAJ_CATA, FIN
  File "C:\work\miniforge3\envs\codeaster-deps\Library\lib\aster\code_aster\Commands.py", line 29, in <module>
    from .Utilities.rc import rc
  File "C:\work\miniforge3\envs\codeaster-deps\Library\lib\aster\code_aster\Utilities\__init__.py", line 59, in <module>
    from .ExecutionParameter import ExecutionParameter
  File "C:\work\miniforge3\envs\codeaster-deps\Library\lib\aster\code_aster\Utilities\ExecutionParameter.py", line 50, in <module>
    import libaster
ImportError: DLL load failed while importing libaster: The specified module could not be found.
To run manually, use:
. C:\work\miniforge3\envs\codeaster-deps\Library\share\aster\profile.bat
C:\work\Miniforge3\envs\codeaster-deps\python.exe C:\work\code\code-aster-src\build\std\debug\catalo\fort.1 --memory=5120
```

Now, this is caused by the current clash between llvm-openmp (needed by intel-fortran-rt) and 
intel-openmp (required by mkl) See https://github.com/conda-forge/conda-forge.github.io/issues/1597.

To resolve this issue, the following command can be used:

```
mamba remove llvm-openmp intel-openmp --force
```

Then, reinstall only intel-openmp:

```
mamba install intel-openmp
```

Then re-run the compilation. The error message should now be resolved.


