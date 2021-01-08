.. _devguide-lowlevel:


*********************
Low level programming
*********************

=============
Configuration
=============

This section gives some details about the configuration of the source code.
This task is performed by executing ``./waf configure`` in the source tree.


Precompilation variables
------------------------

Precompilation variables are always named in uppercase letters by convention.

To tell that a feature is available, a variable ``HAVE_<feature-name>`` is
to ``1``.
To tell that a feature is disabled, we set ``DISABLE_<feature-name>`` to ``1``.
But **only one is defined**, depending of the default behavior.

Examples: ``HAVE_PETSC`` (feature is installed), ``DISABLE_MPI_CHECK`` (in most
cases it is enabled).

List of supported variables:

- ``DISABLE_MATHLIB_FPE``
- ``DISABLE_MPI_CHECK``
- ``HAVE_BACKTRACE``
- ``HAVE_GETLINE`` (only for Metis interface)
- ``HAVE_HDF5``
- ``HAVE_MED``
- ``HAVE_METIS``
- ``HAVE_MFRONT``
- ``HAVE_MPI``
- ``HAVE_MUMPS``
- ``HAVE_PETSC4PY``
- ``HAVE_PETSC``
- ``HAVE_SCOTCH``
- ``HAVE_TRACEBACKQQ``
