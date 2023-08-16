.. _devguide-pybind:

**********************************
Python wrapping of *DataStructure*
**********************************


Files structures
================

The definition of exposed methods of the *DataStructure* objects is included
into the ``libaster.so``. The main file is
:file:`bibcxx/PythonBindings/LibAster.cxx` that exports all symbols.

By convention, one file is created per derivated *DataStructure* object named
:file:`<DataStructure object>Interface.{{h,cxx}}`.


Exported methods
================

Constructors
------------

*DataStructure* objects must be wrapped by a ``std::shared_ptr``.

The definition of the wrapper must reflect the inheritance of the underlying
C++ DataStructures. In the following example, ``Function`` is derivated
from ``BaseFunction``. It is necessary to pass a Python ``Function``
object where a generic ``BaseFunction`` is expected.

.. note:: The constructors of the underlying *instance* should not be available
    to the final user. That's why the definition should use ``py::no_init``.

To allow the serialization of the *DataStructure* objects using the Python
pickling mechanism, we need a constructor that accepts the *Jeveux* name.
The final user must not call this constructor.

For the creation of the *DataStructure* from Python, we also need a *default*
constructor (with or without any argument).

.. note:: These constructors are defined using factory functions by
    templating ``initFactoryPtr``.


Examples
--------

In this simple example one defines a *default* constructor (without argument)
and the constructor used during unpickling that accepts the *Jeveux* name:

.. code-block:: c++

    void exportFunctionToPython( py::module_ &mod )
    {
        namespace py = pybind11;
        ...

        py::class_< Function, Function::FunctionPtr, BaseFunction > (
            mod, "Function" )
            .def( "__init__", py::make_constructor(
                &initFactoryPtr< Function >) )
            .def( "__init__", py::make_constructor(
                &initFactoryPtr< Function,
                                 std::string >) )
            .def( "setValues", &Function::setValues )
            ...
        ;
    }

A more complex example where the constructor needs another object. So the
*default* constructor needs a pointer on a :class:`~code_aster.Objects.Model`
and another that takes the *Jeveux* name and this pointer:

.. code-block:: c++

    void exportElementaryCharacteristicsToPython( py::module_ &mod )
    {
        using namespace pybind11;

        class_< ElementaryCharacteristics, ElementaryCharacteristics::ElementaryCharacteristicsPtr,
                bases< DataStructure > > ( "ElementaryCharacteristics", no_init )
            .def( "__init__", make_constructor(
                &initFactoryPtr< ElementaryCharacteristics,
                                 ModelPtr >) )
            .def( "__init__", make_constructor(
                &initFactoryPtr< ElementaryCharacteristics,
                                 std::string,
                                 ModelPtr >) )
        ;
    };


Other methods
-------------

See :ref:`devguide-recommendations` for methods with default arguments.


Pickling support
================

See :py:mod:`code_aster.Supervis.Serializer` module for the serialization
management.

- DataStructures usually delegate the serialization to their Python objects.

- Constructors arguments defined by :py:meth:`__getinitargs__` implemented in
  :py:mod:`code_aster.Objects.DataStructure` for most of the classes.

  Example: :py:class:`~code_aster.Objects.ElementaryCharacteristics` defines
  its own arguments.
  :py:meth:`~code_aster.Objects.ElementaryCharacteristics.__getinitargs__`
  returns a tuple with two elements: the *Jeveux* object name and the
  :py:class:`~code_aster.Objects.Model` that are passed to the constructor.

- To restore the internal state of the object, subclasses should defined their
  own :py:meth:`__getstate__` and :py:meth:`__setstate__` methods.
  This is done through a subclass of
  :py:class:`~code_aster.Objects.Serialization.InternalStateBuilder`.

  Example: :py:class:`~code_aster.Objects.AssemblyMatrix` does not take its
  :py:class:`~code_aster.Objects.DOFNumbering` as argument in its constructor.
  So it is stored by the ``save()`` method of its
  :py:class:`~code_aster.Objects.Serialization.InternalStateBuilder` and
  restored by the ``restore()`` method.
