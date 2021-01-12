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

To tell that a feature is available, a variable ``ASTER_HAVE_<feature-name>`` is
to ``1``.
To tell that a feature is disabled, we set ``ASTER_DISABLE_<feature-name>`` to ``1``.
But **only one is used**, depending of the default behavior.

If ``ASTER_HAVE_xxx`` is defined (don't check the value, use ``#ifdef``, not
``#if ASTER_HAVE_xxx == 1``), the feature *xxx* is available.
The exception is for ``ASTER_HAVE_PETSC4PY`` because Cython requires that the variable
is always defined.

Examples: ``ASTER_HAVE_PETSC`` (feature is installed), ``ASTER_DISABLE_MPI_CHECK`` (in most
cases it is enabled).

List of supported variables:

- ``ASTER_HAVE_SUPPORT_FPE``
- ``ASTER_DISABLE_MPI_CHECK``
- ``ASTER_HAVE_BACKTRACE``
- ``ASTER_HAVE_GETLINE`` (only for Metis interface)
- ``ASTER_HAVE_HDF5``
- ``ASTER_HAVE_MED``
- ``ASTER_HAVE_METIS``
- ``ASTER_HAVE_MFRONT``
- ``ASTER_HAVE_MPI``
- ``ASTER_HAVE_MUMPS``
- ``ASTER_HAVE_PETSC4PY``
- ``ASTER_HAVE_PETSC``
- ``ASTER_HAVE_SCOTCH``
- ``ASTER_HAVE_TRACEBACKQQ``
