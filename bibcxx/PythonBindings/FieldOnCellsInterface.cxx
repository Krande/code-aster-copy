/**
 * @file FieldOnCellsInterface.cxx
 * @brief Interface python de FieldOnCells
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Code_Aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.
 */

/* person_in_charge: nicolas.sellenet at edf.fr */

#include <boost/python.hpp>

namespace py = boost::python;
#include <PythonBindings/factory.h>

#include "PythonBindings/DataStructureInterface.h"
#include "PythonBindings/FieldOnCellsInterface.h"

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( print_overloads, printMedFile, 1, 2 )

void exportFieldOnCellsToPython() {
    py::class_< FieldOnCellsReal, FieldOnCellsRealPtr, py::bases< DataField > >( "FieldOnCellsReal",
                                                                                 py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< FieldOnCellsReal > ) )
        .def( "__init__", py::make_constructor( &initFactoryPtr< FieldOnCellsReal, std::string > ) )
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< FieldOnCellsReal, ModelPtr, BehaviourPropertyPtr,
                                               std::string, ElementaryCharacteristicsPtr >,
                              py::default_call_policies(),
                              ( py::arg( "model" ), py::arg( "behaviour" ), py::arg( "typcham" ),
                                py::arg( "carael" ) = py::object() ) ) )
        .def( py::init< const FieldOnCellsReal & >() )
        .def( "duplicate", &FieldOnCellsReal::duplicate,R"(
                  Return a duplicated FieldOnCellsReal as a copy
                  Returns:
                  FieldOnCellsReal
                  )",
              ( py::arg( "self" ) ) )
        .def( "exportToSimpleFieldOnCells", &FieldOnCellsReal::exportToSimpleFieldOnCells )
        .def( "getModel", &FieldOnCellsReal::getModel, R"(
                  Return the model associated with the FieldOnCellsReal object
                  Returns:
                  Model: Model Object
                  )",
              ( py::arg( "self" ) )  )
        .def( "getMesh", &FieldOnCellsReal::getMesh,  R"(
                  Return the Mesh associated with the FieldOnCellsReal object
                  Returns:
                  BaseMesh: Mesh object
                  )",
              ( py::arg( "self" ) ) )
        .def( "setDescription", &FieldOnCellsReal::setDescription )
        .def( "setModel", &FieldOnCellsReal::setModel )
        .def( "build", &FieldOnCellsReal::build )
        .def(
            "__getitem__", +[]( const FieldOnCellsReal &v, ASTERINTEGER i ) { return v[i]; } )
        .def(
            "__setitem__",
            +[]( FieldOnCellsReal &v, ASTERINTEGER i, float f ) { return v.operator[]( i ) = f; } )
        .def(
            "__len__", +[]( const FieldOnCellsReal &v ) { return v.size(); } )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self * float() )
        .def( float() * py::self )
        .def( -py::self )
        .def( "setValues", &FieldOnCellsReal::setValues,
              R"(
Set values of the field

Arguments:
    float: value to set
                  )",
              ( py::arg( "self" ), py::arg( "value" ) ) )
        .def( "getValues", &FieldOnCellsReal::getValues,
              py::return_value_policy< py::copy_const_reference >(), R"(
Return a list of values as (x1, y1, z1, x2, y2, z2...)

Returns:
    list[float]: List of values.
                  )",
              ( py::arg( "self" ) ) )
        .def( "size", &FieldOnCellsReal::size, R"(
Return the size of the field

Returns:
    int: number of element in the field
                        )",
              ( py::arg( "self" ) ) )
        .def( "transform", &FieldOnCellsReal::transform< ASTERDOUBLE >, R"(
Apply Function to each value of _ValuesList of the FieldOnCells object.

Arguments:
    func (Function): Callable Python object

Returns:
    FieldOnCellsReal: New FieldOnCells object with the trasformed values
                        )",
              ( py::arg( "self" ), py::arg( "func" ) ) )
        .def( "printMedFile", &FieldOnCellsReal::printMedFile, print_overloads( R"(
Print the field in MED format.

Arguments:
    filename (str): Path to the file to be printed.
    local (bool): Print local values only (relevant for ParallelMesh only, default: *True*)

Returns:
    bool: *True* if succeeds, *False* otherwise.
                        )",
              ( py::arg( "self" ), py::arg( "filename" ), py::arg( "local" ) ) ) )
        .def( "norm", &FieldOnCellsReal::norm< ASTERDOUBLE >, R"(
Return the euclidean norm of the field

Arguments:
    normType (str): "NORM_1", "NORM_2", "NORM_INFINITY"

Returns:
    float: euclidean norm
                        )" )
        .def( "dot", &FieldOnCellsReal::dot< ASTERDOUBLE >, R"(
Return the dot product of two fields

Arguments:
    field (FieldOnNodes): other field

Returns:
    float: dot product
                              )",
              ( py::arg( "self" ), py::arg( "other" ) ) );

    py::class_< FieldOnCellsComplex, FieldOnCellsComplexPtr,
            py::bases< DataField > >( "FieldOnCellsComplex", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< FieldOnCellsComplex >) )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< FieldOnCellsComplex, std::string >) )
        .def(py::init<const FieldOnCellsComplex&>() )
        .def( "duplicate", &FieldOnCellsComplex::duplicate)
        .def( "getModel", &FieldOnCellsComplex::getModel )
        .def( "setDescription", &FieldOnCellsComplex::setDescription )
        .def( "setModel", &FieldOnCellsComplex::setModel )
        .def( "build", &FieldOnCellsComplex::build )
        .def( "getValues", &FieldOnCellsComplex::getValues,
               py::return_value_policy<py::copy_const_reference>(),R"(
                  Return a list of values as (x1, y1, z1, x2, y2, z2...)
                  Returns:
                  list[complex]: List of values.
                  )",
              ( py::arg( "self" ) ))
        .def( "__getitem__",
              +[]( const FieldOnCellsComplex& v, int i ) { return v[i]; } )
        .def( "__setitem__",
              +[]( FieldOnCellsComplex &v, ASTERINTEGER i, ASTERCOMPLEX f )
               { return v.operator[]( i )=f; } )
        .def( "__len__",
              +[]( const FieldOnCellsComplex& v ) { return v.size(); } )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self * float() )
        .def( float() * py::self )
        .def( "size", &FieldOnCellsComplex::size, R"(
Return the size of the field

Returns:
    int: number of element in the field
                        )",
                              ( py::arg( "self" ) ) )
        .def( "transform", &FieldOnCellsComplex::transform<ASTERCOMPLEX>, R"(
Apply Function to each value of _ValuesList of the FieldOnCells object.

Arguments:
    func (Function): Callable Python object

Returns:
    bool: New FieldOnCells object with the trasformed values
                        )",
                              ( py::arg( "self" ), py::arg( "func" ) ))
        .def( "printMedFile", &FieldOnCellsComplex::printMedFile, R"(
Print the field in MED format.

Arguments:
    filename (str): Path to the file to be printed.

Returns:
    bool: *True* if succeeds, *False* otherwise.
                        )",
                  ( py::arg( "self" ), py::arg( "filename" ) )  );


    py::class_< FieldOnCellsLong, FieldOnCellsLongPtr,
            py::bases< DataField > >( "FieldOnCellsLong", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< FieldOnCellsLong >) )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< FieldOnCellsLong, std::string >) )
        .def(py::init<const FieldOnCellsLong&>() )
        .def( "duplicate", &FieldOnCellsLong::duplicate)
        .def( "getModel", &FieldOnCellsLong::getModel )
        .def( "setDescription", &FieldOnCellsLong::setDescription )
        .def( "setModel", &FieldOnCellsLong::setModel )
        .def( "build", &FieldOnCellsLong::build )
        .def( "getValues", &FieldOnCellsLong::getValues,
               py::return_value_policy<py::copy_const_reference>(),R"(
                  Return a list of values as (x1, y1, z1, x2, y2, z2...)
                  Returns:
                  list[int]: List of values.
                  )",
              ( py::arg( "self" ) ))
        .def( "__getitem__",
              +[]( const FieldOnCellsLong& v, int i ) { return v[i]; } )
        .def( "__setitem__",
              +[]( FieldOnCellsLong &v, ASTERINTEGER i, ASTERINTEGER f )
               { return v.operator[]( i )=f; } )
        .def( "__len__",
              +[]( const FieldOnCellsLong& v ) { return v.size(); } )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self * float() )
        .def( float() * py::self )
        .def( "size", &FieldOnCellsLong::size, R"(
Return the size of the field

Returns:
    int: number of element in the field
                        )",
                              ( py::arg( "self" ) ) )
        .def( "printMedFile", &FieldOnCellsLong::printMedFile, R"(
Print the field in MED format.

Arguments:
    filename (str): Path to the file to be printed.

Returns:
    bool: *True* if succeeds, *False* otherwise.
                        )",
                  ( py::arg( "self" ), py::arg( "filename" ) )  );
};
