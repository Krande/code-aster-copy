.. _devguide-debugging:

***********************
Debugging and profiling
***********************


Debugging
=========

This section gives some tips about debugging code_aster.

.. note::

    ``waf_debug install`` turns on ``ASTER_DEBUG_CXX`` flag that prints some
    debugging informations.
    Additional details about objects live-cycle can be printed by defining
    environment variables ``CXXFLAGS=-DASTER_DEBUG_CXX_OBJECTS`` or
    ``CXXFLAGS=-DASTER_DEBUG_CXX_LOW_LEVEL``, or both
    ``CXXFLAGS='-DASTER_DEBUG_CXX_OBJECTS -DASTER_DEBUG_CXX_LOW_LEVEL''``.
    Please ask a guru before adding new blocks under these flags!

    Detailed informations about the command syntax checker can be printed using
    the debug level of the :py:class:`~code_aster.Utilities.logger`.
    Use ``--debug`` option or set ``DEBUG=1`` environment variable.


.. todo::

    It will give an example of debugging the C++/Fortran objects of a
    sequential or a parallel build.

    Another part will give some informations to debug the Python part.


Helper functions
~~~~~~~~~~~~~~~~

- :py:func:`~code_aster.Objects.DataStructure.debugPrint`:
  This method is available on all :py:class:`~code_aster.Objects.DataStructure`
  objects. It prints the content of all its *Jeveux* objects.

- :py:func:`~code_aster.Objects.DataStructure.use_count`:
  This method is a wrapping to ``std::shared_ptr< T >`` ``use_count()``
  member. It is available for only some
  :py:class:`~code_aster.Objects.DataStructure`.


Build of elements failed
~~~~~~~~~~~~~~~~~~~~~~~~

In the case that the installation (``waf_debug install``) failed during building
the elements catalog, the output ends with something like:

.. code-block:: none

    + build the elements catalog elem.1 using installed aster (from cata_ele.ojb)
    stdout: ...
    <sometimes the message is not clear enough>
    ...
    stderr: ...
    To run manually, use:
    ...
    <command line to reproduce the error>

First, you can restart in verbose mode (and with one process to have a synced output):

.. code-block:: sh

    ./waf_debug install -v -j 1

If it is not sufficient, try to reproduce the error in a new terminal (adjust build
and installation directories if needed):

.. code-block:: shell

    cd /tmp
    ulimit -c unlimited
    . ${HOME}/dev/codeaster/install/debug/share/aster/profile.sh
    python3 ${HOME}/dev/codeaster/src/build/mpidebug/debug/catalo/fort.1 --memory=4096

It may create a core file. In this case, try:

.. code-block:: shell

    gdb $(which python3) core
    (gdb) where
    ...

If it does not give interesting informations, try to import each modules of
code_aster as done in :file:`code_aster/__init__.py`:

.. code-block:: python

    python3
    >>> import aster
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    ImportError: .../dev/codeaster/install/mpi/lib64/aster/libbibfor.so: undefined symbol: scotchfdgraphcorderinit_

This is an example of error caused by a missing external library.

Another frequent error:

.. code-block:: python

    python3
    >>> import aster
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    ImportError: .../dev/codeaster/install/debug/lib64/aster/libbibcxx.so: undefined symbol: _Z7nmdocr_PKcS0_Pcjj

Here, ``nmdocr`` is a Fortran subroutine, called from C++. Its prototype must be
enclosed by ``extern "C" { ... }``.


Profiling
=========

The well known tool ``gprof`` is a very good and simple choice to profile an
executable but it does not work to profile a shared library.
And code_aster is a Python module built as a shared library.

.. note::

    Profiling code_aster using `gperftools`_ has been tested but the analysis
    of the results was difficult.

    More tools have to be evaluated.


.. _gperftools: https://github.com/gperftools/gperftools
