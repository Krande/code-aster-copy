.. _devguide-recommendations:

***********************
Recommendations / howto
***********************


Definition ``.h`` vs Implementation ``.cxx``
============================================

The header files ``.h`` should only contain definitions. Exceptions are allowed
for *properties* that directly return the content of an attribute
(example: ``DataStructure.getName()``).
Others functions should be implemented in the ``.cxx`` file
(example: ``DataStructure.addDependency()``).

It avoids to rebuild everything for a small change.


Default arguments and overloaded methods in Python Interface
============================================================

When a C++ method has default arguments, the default values must also be explicitly
defined in the Python interface (because they are not part of the function’s type
information). Do not use overloading just for default arguments.

Overloaded methods have to be described in the Python interface using
``py::overload_cast< types... >( method )``.
Define overloaded methods in the Python interface **only** if the types of the arguments
differ.

See as example :meth:`Mesh.getNodes` and its interface in
:file:`bibcxx/PythonBindings/MeshInterface.cxx`.

.. code-block:: c++

    #include "aster_pybind.h";

    void exportMeshToPython( py::module_ &mod )
    {
        py::class_< Mesh, Mesh::MeshPtr, BaseMesh >( mod, "Mesh" )
            .def( "getNodes",
                py::overload_cast< const std::string, const bool, const bool >( &Mesh::getNodes,
                                                                                py::const_ ),
                R"(
                ... docstring ...

                Arguments:
                    group_name (str): Description
                    ...

                Returns:
                    type: Description.
                )",
                py::arg( "group_name" ) = "", py::arg( "localNumbering" ) = true,
                py::arg( "same_rank" ) = true )
            .def( "getNodes",
                py::overload_cast< const bool, const bool >( &Mesh::getNodes, py::const_ ),
                py::arg( "localNumbering" ), py::arg( "same_rank" ) = true )

Note:
    - The docstring should not be repeated, or only if the arguments are very different
      between overloaded methods because all docstrings will be exported in Sphinx pages.
    - Default values must be consistent with those of the class header.
    - ``const`` methods must use the additional argument ``py::const_``.
    - C++ methods should not create ``PyObject`` but ``std::map``, ``std::vector``...


How to return different types from a unique Python interface
============================================================

The situation is as follows:
On the C++ side, we have 2 member functions of the class ``ModeResult``
with different names and which return different types:

.. code-block:: c++

    AssemblyMatrixDisplacementRealPtr getDisplacementRealStiffnessMatrix();
    AssemblyMatrixTemperatureRealPtr getTemperatureRealStiffnessMatrix();


We want to interface these 2 functions by only one ``getStiffnessMatrix`` on the
Python side but which returns the 2 different types according to the case.

To do this, we use the ``std::variant`` type to store both types returned by the 2 functions.

.. code-block:: c++

    using MatrixVariant = std::variant< AssemblyMatrixDisplacementRealPtr,
                                        AssemblyMatrixTemperatureRealPtr >;


Then, we have to write a ``getStiffnessMatrix`` function that returns one or the
other type depending on the case but storing it in a ``std::variant``.
This function will become a member function of the python class
:py:class:`code_aster.Objects.ModeResult` so the C++ function
must take as argument a ``ModeResultPtr``.

.. code-block:: c++

    MatrixVariant getStiffnessMatrix( const ModeResultPtr self )
    {
        auto mat1 = self->getDisplacementRealStiffnessMatrix();
        if( mat1 != nullptr )
            return MatrixVariant( mat1 );
        auto mat2 = self->getTemperatureRealStiffnessMatrix();
        return MatrixVariant( mat2 );
    };


In the pybind11 interface of the ``ModeResult`` class, we must add the function:

.. code-block:: c++

    .def( "getStiffnessMatrix", &getStiffnessMatrix )


NB: In the real life, ``getStiffnessMatrix`` is a template function.


Macro-Commands
==============

Legacy Macro-commands do not work as is.

#. There is no need to define an executor manually.
   Default :class:`~code_aster.Supervis.ExecuteCommand.ExecuteMacro` is just
   adapted by :mod:`code_aster.Commands.operator` using the right catalog
   description.

#. The body of the macro-command, the ``ops()`` function, is automatically
   called by the :meth:`~code_aster.Supervis.ExecuteCommand.ExecuteMacro.run`
   factory.

#. Results of macro-commands are created directly by the ``ops()`` function
   (called by ``exec_()``). ``create_result()`` method does nothing else
   registering the additional results (declared with ``CO()``).

#. The ``ops()`` function must now returns the result object it creates.


For user Macro-commands or those from *Contrib* directory, an executor must be
manually added (since the catalog description can not be imported from the
official ones). A convenient function allows to easily define this executor:

.. code-block:: python

    from code_aster.Supervis.ExecuteCommand import UserMacro
    MA_MACRO = UserMacro("MA_MACRO", MA_MACRO_cata, ma_macro_ops)


Required changes
----------------

- The ``ops()`` function returned an exit code as integer.

  Now, it must return the created result object, or *None* if there is not.

- In code_aster legacy the keywords arguments passed to ``ops()`` contained
  all existing keywords, eventually with *None* value.

  Now, only the user keywords + the default keywords are passed.
  So, only compulsory keywords and those having a default value can be arguments
  of the ``ops()`` function.
  If needed, these arguments may be wrapped by ``_F()`` that provides a ``[]``
  operator that returns *None* if a keyword does not exist.

  Example:

  .. code-block:: python

        def my_macro_ops(INFO, **kwargs):
            """..."""
            kwargs = _F(kwargs)
            para = kwargs['NOM_PARA']  # no failure even if the keyword does not exist

- Tests on DataStructures types must be changed.
  For example:

  Replace ``AsType(obj) is fonction_sdaster``, ``type(obj) is fonction_sdaster``
  or ``isinstance(obj, fonction_sdaster)``

  by ``obj.getType() == "FONCTION"``

- Object ``MCLIST`` does not exist anymore. List of factor keywords is just a
  *list* or a *tuple*.

  Just use :func:`~code_aster.Utilities.force_list` to ensure to have a list
  even if the user passed only one occurrence.

- ``.List_F()`` does not exist anymore.

  Replace ``POUTRE.List_F()`` by ``force_list(POUTRE)``.

  Temporarly one can use ``POUTRE = ListFact(POUTRE)`` not to change the code
  and let ``POUTRE.List_F()`` with a dummy ``.List_F()`` function that does nothing.

- Usage of logical units: See :mod:`code_aster.Helpers.LogicalUnit`.

- Additional results (**CO()** objects):

  They must be registered with
  :meth:`~code_aster.Supervis.ExecuteCommand.ExecuteMacro.register_result`.
  It replaces *DeclareOut()* but must be called **after** the result creation.

  .. code-block:: diff

        -          self.DeclareOut('num', numeddl)
        +          # self.DeclareOut('num', numeddl)
                   num = NUME_DDL(MATR_RIGI=_a, INFO=info)
        +          self.register_result(num, numeddl)

  In the legacy version some testcases sometimes define ``OBJ = CO('NAME')`` and
  then pass either ``NAME`` or ``OBJ`` to children commands.
  Now using the legacy mode of macro-commands that publishes ``NAME`` in the parent
  context ``OBJ`` can not be passed to children commands. It will not have the
  expected type (it stays a ``CO`` object and not becomes a ``Table`` or
  ``Mesh``!).

  When the new mode will be enabled one will just use ``result.NAME`` without
  ambiguity.


Parallel specific DataStructures
================================

Q: How to pass a :py:class:`code_aster.Objects.ParallelMesh` to a command?

A: The solution is in "a :py:class:`code_aster.Objects.ParallelMesh` is a :py:class:`code_aster.Objects.Mesh`". It is just necessary to declare a
DataStructure is the Python command description (*catalog*) that matches the
same type.
Example: :py:meth:`code_aster.Objects.ParallelMesh.getType()`
returns ``MAILLAGE_P``, so one defines:

.. code-block:: python

    class maillage_p(maillage_sdaster):
        pass


Common errors
=============

- The compilation works but ``waf_debug install`` ends with
  ``stderr: Segmentation fault`` during the compilation of elements catalogs.

  **Explanation**: It may be an error in a Python function called from a C or
  Fortran function.
  Check it by manually importing the module in a Python interpreter:

  .. code-block:: sh

      $ cp ../src/build/mpidebug/catalo/cata_ele.ojb fort.4
      $ python
      >>> from code_aster import CA
      >>> CA.init(CATALOGUE={"FICHIER": "CATAELEM", "UNITE": 4})
      >>> CA.MAJ_CATA(ELEMENT={})
      >>> exit()

  An undefined symbol in an underlying library, for example ``libaster.so``,
  may also cause ``stderr: Segmentation fault``.
  Try to import the libraries one by one.

- The compilation works but the first execution fails with
  ``ImportError: generic_type: type "PythonBool" is already registered!``.

  **Solution**: You just adds a pybind11 binding for a new class?
  You probably forgot the ``shared::ptr`` ancestor (``xxPtr``).

  .. code-block:: c++

    py::class_< Picklable, PicklablePtr >( mod, "Picklable" )
