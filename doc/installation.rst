.. _building:


##################################
Building and installing code_aster
##################################


*******************
Getting the sources
*******************

The code_aster sources are under `Git <https://git-scm.com/>`_
version control system.

Use these commands to check out the latest project source code from the
`GitLab repository <https://gitlab.com/codeaster/src>`_::

    $ mkdir -p $HOME/dev/codeaster
    $ cd $HOME/dev/codeaster
    $ git clone https://gitlab.com/codeaster/src.git

In this manual, we assume that the working directory corresponds to the code_aster
source folder (``cd $HOME/dev/codeaster/src``).


********************
Compiling code_aster
********************

The building and installation instructions can be found on documentation repository
`Installation and Development
<https://gitlab.com/codeaster-opensource-documentation/opensource-installation-development>`_.


.. _testing:

*************************
Testing your installation
*************************

The :py:mod:`run_aster` package provides a convenient script that wraps ``ctest`` to
execute a list of testcases.

Examples:

- Check the sequential testcases of the ``submit`` list (~2000 tests) using a parallel build:

  .. code-block:: sh

    ../install/mpi/bin/run_ctest --resutest=../resutest -L submit -L sequential

- Check a selection of testcases:

  .. code-block:: sh

    cat << EOF > /tmp/some_tests
    ssll112a
    ssnl125b
    zzzz144a
    zzzz388a
    EOF
    ../install/mpi/bin/run_ctest --resutest=../resutest --testlist=/tmp/some_tests

If you are using a sequential build of code_aster and that the testlist contains
parallel testcases, they are skipped.


**********************
Advanced configuration
**********************

- ``addmem`` parameter

  code_aster tries not to use more memory than the value passed to the
  ``--memory`` option. This memory is allocated for the objects of the study,
  but also for the executable itself and the loaded shared libraries.

  The ``addmem`` value (in MB) represents the memory used to load the executable
  and its libraries. A good evaluation is obtained by just running a minimal
  ``DEBUT()``/``FIN()`` execution and seing the value printed at
  *"MAXIMUM DE MEMOIRE UTILISEE PAR LE PROCESSUS"*.

  The value passed to the ``--memory`` option is equal to ``addmem`` +
  *the value requested by the user* in the graphical user interface.

  .. note:: This value can be more than 2000 MB using Intel MPI for example.

  See :py:mod:`run_aster.config` for more information.
