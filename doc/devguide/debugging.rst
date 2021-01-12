.. _devguide-debugging:

***********************
Debugging and profiling
***********************


Debugging
=========

This section gives some tips about debugging code_aster.

.. note::

    ``waf install_debug`` turns on ``ASTER_DEBUG_CXX`` flag that prints some
    informations about objects live-cycle.
    Please ask a guru before adding new blocks under this flag!

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
  This method is a wrapping to ``boost::shared_ptr< T >`` ``use_count()``
  member. It is available for only some
  :py:class:`~code_aster.Objects.DataStructure`.


Build of elements failed
~~~~~~~~~~~~~~~~~~~~~~~~

In the case that the installation (``waf install_debug``) failed during building
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


Try to reproduce the error in a new terminal (adjust ``std`` to ``mpi`` if needed):

.. code-block:: shell

    cd /tmp
    ulimit -c unlimited
    . ${HOME}/dev/codeaster/install/std/share/aster/profile.sh
    python3 ${HOME}/dev/codeaster/src/build/std/debug/catalo/fort.1 --memory=4096

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

(This is an example of error!)



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
