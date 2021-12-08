/**
 * @file FieldOnNodesInterface.cxx
 * @brief Interface python de FieldOnNodes
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

#include "DataFields/MeshCoordinatesField.h"
#include "PythonBindings/DataStructureInterface.h"
#include "PythonBindings/FieldOnNodesInterface.h"

void exportFieldOnNodesToPython() {

    py::class_< FieldOnNodesReal, FieldOnNodesRealPtr,
                py::bases< DataField > >( "FieldOnNodesReal", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< FieldOnNodesReal >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< FieldOnNodesReal, std::string >))
        .def(py::init<const FieldOnNodesReal&>())
        .def( "duplicate", &FieldOnNodesReal::duplicate)
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< FieldOnNodesReal,
                                                    BaseDOFNumberingPtr>))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< FieldOnNodesReal,
                                                    BaseDOFNumberingPtr, JeveuxMemory >))
        .def( "exportToSimpleFieldOnNodes",
              &FieldOnNodesReal::exportToSimpleFieldOnNodes )
        .def( "getMesh", &FieldOnNodesReal::getMesh )
        .def( "__getitem__",
              +[]( const FieldOnNodesReal &v, int i ) { return v.operator[]( i ); } )
        .def( "__setitem__",
              +[]( FieldOnNodesReal &v, int i, float f ) { return v.operator[]( i )=f; } )
        .def( "__add__",
              +[]( FieldOnNodesReal &v1, FieldOnNodesReal &v2 ) { return (v1 + v2); } )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( py::self * float() )
        .def( float() * py::self )
        .def( py::self *= float() )
        .def( - py::self)
        .def( "printMedFile", &FieldOnNodesReal::printMedFile )
        .def( "setDOFNumbering", &FieldOnNodesReal::setDOFNumbering )
        .def( "setMesh", &FieldOnNodesReal::setMesh )
        .def( "setDescription", &FieldOnNodesReal::setDescription )
        .def( "build", &FieldOnNodesReal::build )
        .def( "getDOFNumbering", &FieldOnNodesReal::getDOFNumbering )
        .def( "getMesh", &FieldOnNodesReal::getMesh )
        .def( "getDescription", &FieldOnNodesReal::getDescription )
        .def( "updateValuePointers", &FieldOnNodesReal::updateValuePointers )
        .def( "norm", &FieldOnNodesReal::norm< ASTERDOUBLE >,
              R"(
Return the euclidean norm of the field

Argument:
    normType: "NORM_1", "NORM_2", "NORM_INFINITY"

Returns:
    double: euclidean norm
        )",
              ( py::arg( "self" ) ) )
        .def( "dot", &FieldOnNodesReal::dot< ASTERDOUBLE >,
              R"(
Return the dot product of two fields

Argument:
    FieldOnNodes: field

Returns:
    double: dot produc
        )",
              ( py::arg( "self"), py::arg( "field") ) )
        .def( "size", &FieldOnNodesReal::size,
               R"(
Return the size of the field

Return:
    int: number of element in the field
        )",
              ( py::arg( "self" ) ) )
        .def( "setValues", &FieldOnNodesReal::setValues,
               R"(
Set values of the field

Argument:
    float: value to set
        )",
              ( py::arg( "self" ), py::arg( "value" ) ) )
        .def( "getValues", &FieldOnNodesReal::getValues,
        py::return_value_policy<py::copy_const_reference>(), R"(
Return a list of values as (x1, y1, z1, x2, y2, z2...)

Returns:
    list[float]: List of values.
        )",
              ( py::arg( "self" ) ) );
              
    py::class_< FieldOnNodesComplex, FieldOnNodesComplexPtr,
                py::bases< DataField > >( "FieldOnNodesComplex", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< FieldOnNodesComplex >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< FieldOnNodesComplex, std::string >))
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< FieldOnNodesComplex, BaseDOFNumberingPtr > ) )
        .def( "exportToSimpleFieldOnNodes",
              &FieldOnNodesComplex::exportToSimpleFieldOnNodes )
        .def( "getMesh", &FieldOnNodesComplex::getMesh )
        .def( "__getitem__",
              +[]( const FieldOnNodesComplex &v, int i ) { return v.operator[]( i ); } )
        .def(
            "__setitem__",
            +[]( FieldOnNodesComplex &v, int i, ASTERCOMPLEX f ) { return v.operator[]( i ) = f; } )
        .def( "printMedFile", &FieldOnNodesComplex::printMedFile )
        .def( "setDOFNumbering", &FieldOnNodesComplex::setDOFNumbering )
        .def( "setMesh", &FieldOnNodesComplex::setMesh )
        .def( "setDescription", &FieldOnNodesComplex::setDescription )
        .def( "build", &FieldOnNodesComplex::build )
        .def( "getDOFNumbering", &FieldOnNodesComplex::getDOFNumbering )
        .def( "getMesh", &FieldOnNodesComplex::getMesh )
        .def( "getDescription", &FieldOnNodesComplex::getDescription )
        .def( "getValues", &FieldOnNodesComplex::getValues,
              py::return_value_policy< py::copy_const_reference >(),R"(
Return a list of complex values as [x11, x21, ..., xm1, x12, x22, ..., xm2...] 
(m is the total number of componenets)

Returns:
    list[complex]: List of values.
        )",
              ( py::arg( "self" ) ) )
        .def( "setValues", &FieldOnNodesComplex::setValues, R"(
Set values of the field

Argument:
    complex: value to set
        )",
              ( py::arg( "self" ), py::arg( "value" ) ) )
        .def( "updateValuePointers", &FieldOnNodesComplex::updateValuePointers );
};
