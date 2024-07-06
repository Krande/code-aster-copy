# Overflow Error


```
# ------------------------------------------------------------------------------
Command line #1:
    python .\comp008e.comm.changed.py --test --last --memory 512.0 --tpmax 1500 --numthreads 1 >> fort.6 2>&1
Traceback (most recent call last):
  File "<frozen runpy>", line 198, in _run_module_as_main
  File "<frozen runpy>", line 88, in _run_code
  File "C:\Work\miniforge3\envs\codeaster-deps\Library\lib\aster\run_aster\run_aster_main.py", line 506, in <module>
    sys.exit(main())
             ^^^^^^
  File "C:\Work\miniforge3\envs\codeaster-deps\Library\lib\aster\run_aster\run_aster_main.py", line 492, in main
    status = calc.execute(wrkdir)
             ^^^^^^^^^^^^^^^^^^^^
  File "C:\Work\miniforge3\envs\codeaster-deps\Library\lib\aster\run_aster\run.py", line 140, in execute
    status = self._execute()
             ^^^^^^^^^^^^^^^
  File "C:\Work\miniforge3\envs\codeaster-deps\Library\lib\aster\run_aster\run.py", line 156, in _execute
    status = self.execute_study()
             ^^^^^^^^^^^^^^^^^^^^
  File "C:\Work\miniforge3\envs\codeaster-deps\Library\lib\aster\run_aster\run.py", line 196, in execute_study
    status.update(self._exec_one(comm, timeout - status.times[-1]))
                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Work\miniforge3\envs\codeaster-deps\Library\lib\aster\run_aster\run.py", line 218, in _exec_one
    exitcode = run_command(cmd, exitcode_file=EXITCODE_FILE)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Work\miniforge3\envs\codeaster-deps\Library\lib\aster\run_aster\utils.py", line 155, in run_command
    iret = waitstatus_to_exitcode(iret)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
OverflowError: can't convert negative int to unsigned
```