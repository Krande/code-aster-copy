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

    py::class_< FieldOnNodesRealClass, FieldOnNodesRealPtr,
                py::bases< DataFieldClass > >( "FieldOnNodesReal", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< FieldOnNodesRealClass >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< FieldOnNodesRealClass, std::string >))
        .def(py::init<const FieldOnNodesRealClass&>())
        .def( "duplicate", &FieldOnNodesRealClass::duplicate)
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< FieldOnNodesRealClass,
                                                    BaseDOFNumberingPtr>))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< FieldOnNodesRealClass,
                                                    BaseDOFNumberingPtr, JeveuxMemory >))
        .def( "exportToSimpleFieldOnNodes",
              &FieldOnNodesRealClass::exportToSimpleFieldOnNodes )
        .def( "getMesh", &FieldOnNodesRealClass::getMesh )
        .def( "__getitem__",
              +[]( const FieldOnNodesRealClass &v, int i ) { return v.operator[]( i ); } )
        .def( "__setitem__",
              +[]( FieldOnNodesRealClass &v, int i, float f ) { return v.operator[]( i )=f; } )
        .def( "__add__",
              +[]( FieldOnNodesRealClass &v1, FieldOnNodesRealClass &v2 ) { return (v1 + v2); } )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( py::self * float() )
        .def( float() * py::self )
        .def( py::self *= float() )
        .def( - py::self)
        .def( "printMedFile", &FieldOnNodesRealClass::printMedFile )
        .def( "setDOFNumbering", &FieldOnNodesRealClass::setDOFNumbering )
        .def( "setMesh", &FieldOnNodesRealClass::setMesh )
        .def( "setDescription", &FieldOnNodesRealClass::setDescription )
        .def( "update", &FieldOnNodesRealClass::update )
        .def( "getDOFNumbering", &FieldOnNodesRealClass::getDOFNumbering )
        .def( "getMesh", &FieldOnNodesRealClass::getMesh )
        .def( "getDescription", &FieldOnNodesRealClass::getDescription )
        .def( "updateValuePointers", &FieldOnNodesRealClass::updateValuePointers )
        .def( "norm", &FieldOnNodesRealClass::norm,
               R"(
Return the euclidean norm of the field

Argument:
    normType: "NORM_1", "NORM_2", "NORM_INFINITY"

Returns:
    double: euclidean norm
        )",
              ( py::arg( "self" ) ) )
        .def( "dot", &FieldOnNodesRealClass::dot,
               R"(
Return the dot product of two fields

Argument:
    FieldOnNodes: field

Returns:
    double: dot produc
        )",
              ( py::arg( "self"), py::arg( "field") ) )
        .def( "size", &FieldOnNodesRealClass::size,
               R"(
Return the size of the field

Return:
    int: number of element in the field
        )",
              ( py::arg( "self" ) ) )
        .def( "setValues", &FieldOnNodesRealClass::setValues,
               R"(
Set values of the field

Argument:
    float: value to set
        )",
              ( py::arg( "self" ), py::arg( "value" ) ) )
        .def( "getValues", &FieldOnNodesRealClass::getValues,
        py::return_value_policy<py::copy_const_reference>(), R"(
Return a list of values as (x1, y1, z1, x2, y2, z2...)

Returns:
    list[float]: List of values.
        )",
              ( py::arg( "self" ) )   );
    py::class_< FieldOnNodesComplexClass, FieldOnNodesComplexPtr,
                py::bases< DataFieldClass > >( "FieldOnNodesComplex", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< FieldOnNodesComplexClass >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< FieldOnNodesComplexClass, std::string >))
        .def( "exportToSimpleFieldOnNodes",
              &FieldOnNodesComplexClass::exportToSimpleFieldOnNodes )
        .def( "getMesh", &FieldOnNodesComplexClass::getMesh )
        .def( "__getitem__",
              +[]( const FieldOnNodesComplexClass &v, int i ) { return v.operator[]( i ); } )
        .def( "printMedFile", &FieldOnNodesComplexClass::printMedFile )
        .def( "setDOFNumbering", &FieldOnNodesComplexClass::setDOFNumbering )
        .def( "setMesh", &FieldOnNodesComplexClass::setMesh )
        .def( "setDescription", &FieldOnNodesComplexClass::setDescription )
        .def( "update", &FieldOnNodesComplexClass::update )
        .def( "getDOFNumbering", &FieldOnNodesComplexClass::getDOFNumbering )
        .def( "getMesh", &FieldOnNodesComplexClass::getMesh )
        .def( "getDescription", &FieldOnNodesComplexClass::getDescription )
        .def( "update", &FieldOnNodesComplexClass::update )
        .def( "updateValuePointers", &FieldOnNodesComplexClass::updateValuePointers );
};
