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
To tell that a feature is disabled, we set ``ASTER_HAVE_<feature-name>`` is undefined.

Examples: ``ASTER_HAVE_PETSC`` (feature is installed).

If ``ASTER_HAVE_xxx`` is defined (use ``#ifdef``, not
``#if ASTER_HAVE_xxx == 1``), the feature *xxx* is available.

Include guards to avoid multiple inclusion must be named ``NAME_OF_HEADER_FILE_H_``.
Example: ``MATERIAL_H_`` for :file:`Material.h`.

Variables about the availability of prerequisites:

- ``ASTER_HAVE_HDF5``
- ``ASTER_HAVE_MED``
- ``ASTER_HAVE_METIS``
- ``ASTER_HAVE_MFRONT``
- ``ASTER_HAVE_MUMPS``
- ``ASTER_HAVE_PETSC``
- ``ASTER_HAVE_PETSC4PY``
- ``ASTER_HAVE_SCOTCH``

Type of compiler and mathematical libraries:

- ``ASTER_HAVE_INTEL_IFORT``
- ``ASTER_HAVE_MPI``
- ``ASTER_HAVE_OPENMP``
- ``ASTER_HAVE_MKL``
- ``ASTER_HAVE_OPENBLAS``

Use in C/Fortran interfaces:

- ``ASTER_ADD_STDCALL``
- ``ASTER_NO_UNDERSCORE``
- ``ASTER_STRINGIFY_USE_OPERATOR``
- ``ASTER_STRINGIFY_USE_QUOTES``
- ``ASTER_STRLEN_AT_END``

Platform types:

- ``ASTER_PLATFORM_DARWIN``
- ``ASTER_PLATFORM_DARWIN64``
- ``ASTER_PLATFORM_FREEBSD``
- ``ASTER_PLATFORM_FREEBSD64``
- ``ASTER_PLATFORM_LINUX``
- ``ASTER_PLATFORM_LINUX64``
- ``ASTER_PLATFORM_POSIX``
- ``ASTER_PLATFORM_WINDOWS``

Datatypes and adjustements to prerequisites datatypes:

- ``ASTER_HAVE_64_BITS``: Enabled on 64 bits OS, use integer8 as default
  integers.
- ``ASTER_LOGICAL_SIZE``
- ``ASTER_MED_SAME_INT_IDT``: *True* if code_aster and Med use same integers
  types (not conversion needed in arguments passing).
- ``ASTER_MED_VERSION_MAJOR``
- ``ASTER_MED_VERSION_MINOR``
- ``ASTER_MULT_FRONT_BLOCK_SIZE``
- ``ASTER_MUMPS_VERSION``
- ``ASTER_MUMPS_CONSORTIUM``
- ``ASTER_MUMPS_REDUCMPI``
- ``ASTER_PETSC_64BIT_INDICES``: Enabled if PETSc uses integer8 for *PetscInt*.
- ``ASTER_PETSC_HAVE_HYPRE``
- ``ASTER_PETSC_HAVE_ML``
- ``ASTER_PETSC_HAVE_MUMPS``
- ``ASTER_PETSC_HAVE_SUPERLU``
- ``ASTER_PETSC_INT_SIZE``: Size of *PetscInt*.

Advanced features:

- ``ASTER_ENABLE_MPI_CHECK``: Enables checking before MPI communications.
- ``ASTER_ENABLE_PROC_STATUS``: Enabled if ``/proc/status`` is available.
- ``ASTER_HAVE_SUPPORT_FPE``: Enabled if FPE can be catched.
- ``ASTER_HAVE_BACKTRACE``: Enabled if GNU traceback function is available.
- ``ASTER_HAVE_TRACEBACKQQ``: Enabled if Intel traceback function is available.
- ``ASTER_IGNORE_WITH_ASLINT``: Locally enabled to skip a source block during
  aslint checkings.
- ``ASTER_TEST_STRICT``: Used to force ``TOLE_MACHINE`` to 1.e-6.
- ``ASTER_WITHOUT_PYMOD``: Used to create a *main* entry.


To print more debugging informations:

- ``ASTER_DEBUG_ALL``
- ``ASTER_DEBUG_ALLOCATE``
- ``ASTER_DEBUG_ASSERT``
- ``ASTER_DEBUG_CXX``
- ``ASTER_DEBUG_CXX_LOW_LEVEL``
- ``ASTER_DEBUG_CXX_OBJECTS``
- ``ASTER_DEBUG_DLL``
- ``ASTER_DEBUG_EXCEPT``
- ``ASTER_DEBUG_FONCTIONS``
- ``ASTER_DEBUG_IODR``
- ``ASTER_DEBUG_LOC``
- ``ASTER_DEBUG_MED``
- ``ASTER_DEBUG_MPI``
- ``ASTER_DEBUG_MPICOM``
- ``ASTER_DEBUG_UNITTEST_MPI``: If enabled, a unittest is run just after
  *MPI_Init*.

System variables used in code_aster source code:

- ``OPEN_MPI``: Enabled if OpenMPI implementation is used.
- ``HAVE_GETLINE``: Used in Metis interface.
- ``NDEBUG``: Standard variable for non-debug build, used to skip code for
  released installations.
- ``NPY_API_VERSION``
