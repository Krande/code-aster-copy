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

- ``ASTER_ADD_STDCALL``

- ``ASTER_DISABLE_MPI_CHECK``
- ``ASTER_ENABLE_PROC_STATUS``
- ``ASTER_HAVE_64_BITS``
- ``ASTER_HAVE_BACKTRACE``
- ``ASTER_HAVE_GETLINE`` (only for Metis interface)
- ``ASTER_HAVE_HDF5``
- ``ASTER_HAVE_INTEL_IFORT``
- ``ASTER_HAVE_MED``
- ``ASTER_HAVE_METIS``
- ``ASTER_HAVE_MFRONT``
- ``ASTER_HAVE_MKL``
- ``ASTER_HAVE_MPI``
- ``ASTER_HAVE_MUMPS``
- ``ASTER_HAVE_OPENBLAS``
- ``ASTER_HAVE_OPENMP``
- ``ASTER_HAVE_PETSC``
- ``ASTER_HAVE_PETSC4PY``
- ``ASTER_HAVE_SCOTCH``
- ``ASTER_HAVE_SUPPORT_FPE``
- ``ASTER_HAVE_TRACEBACKQQ``
- ``ASTER_IGNORE_WITH_ASLINT``
- ``ASTER_LOGICAL_SIZE``
- ``ASTER_MED_SAME_INT_IDT``
- ``ASTER_MED_VERSION_MAJOR``
- ``ASTER_MED_VERSION_MINOR``
- ``ASTER_MULT_FRONT_BLOCK_SIZE``
- ``ASTER_MUMPS_VERSION``
- ``ASTER_NO_UNDERSCORE``
- ``ASTER_PETSC_64BIT_INDICES``
- ``ASTER_PETSC_HAVE_HYPRE``
- ``ASTER_PETSC_HAVE_ML``
- ``ASTER_PETSC_HAVE_MUMPS``
- ``ASTER_PETSC_HAVE_SUPERLU``
- ``ASTER_PETSC_INT_SIZE``
- ``ASTER_PLATFORM_DARWIN``
- ``ASTER_PLATFORM_DARWIN64``
- ``ASTER_PLATFORM_FREEBSD``
- ``ASTER_PLATFORM_FREEBSD64``
- ``ASTER_PLATFORM_LINUX``
- ``ASTER_PLATFORM_LINUX64``
- ``ASTER_PLATFORM_POSIX``
- ``ASTER_PLATFORM_SOLARIS``
- ``ASTER_PLATFORM_SOLARIS64``
- ``ASTER_PLATFORM_WINDOWS``
- ``ASTER_STRINGIFY_USE_OPERATOR``
- ``ASTER_STRINGIFY_USE_QUOTES``
- ``ASTER_STRLEN_AT_END``
- ``ASTER_TEST_STRICT``
- ``ASTER_UNITTEST_MPI``
- ``ASTER_WITHOUT_PYMOD``


To print more debugging informations:

- ``ASTER_DEBUG_ALL``
- ``ASTER_DEBUG_ALLOCATE``
- ``ASTER_DEBUG_ASSERT``
- ``ASTER_DEBUG_CXX``
- ``ASTER_DEBUG_DLL``
- ``ASTER_DEBUG_EXCEPT``
- ``ASTER_DEBUG_FONCTIONS``
- ``ASTER_DEBUG_IODR``
- ``ASTER_DEBUG_LOC``
- ``ASTER_DEBUG_MED``
- ``ASTER_DEBUG_MPI``
- ``ASTER_DEBUG_MPICOM``
