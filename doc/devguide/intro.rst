.. _devguide-general:

#################
General concepts
#################

All C, C++ and fortran source files are compiled and embedded in :file:`libaster.so`.
This shared library also defines Python modules: ``aster``, ``aster_core``,
``aster_fonctions``, ``med_aster`` and ``libaster`` which defines Python bindings
to the C++ objects (using `pybind11 <https://pybind11.readthedocs.io>`_).
The C++ objects may be extended with pure Python methods.

The import of the ``libaster`` module initializes all other Python modules.
The user does not directly import ``libaster`` but the ``code_aster`` package.

The initialization of the memory manager (*Jeveux*) should be explicit done.
The closure of the memory manager is automatically done when there is no more
reference to ``libaster``.

Some Python objects are created during the initialization and are callable
by C/C++ functions (see :py:mod:`code_aster.Supervis` for details).

.. note::
    Please read carefully the :ref:`devguide-recommendations`.
