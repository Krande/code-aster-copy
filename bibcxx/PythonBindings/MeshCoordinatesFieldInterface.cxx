/**
 * @file exportMeshCoordinatesFieldToPython.cxx
 * @brief Interface python de MeshCoordinates
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

#include "PythonBindings/MeshCoordinatesFieldInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/DataStructureInterface.h"
#include "PythonBindings/FieldOnNodesInterface.h"

void exportMeshCoordinatesFieldToPython( py::module_ &mod ) {

    py::class_< MeshCoordinatesField, MeshCoordinatesFieldPtr, DataStructure >(
        mod, "MeshCoordinatesField" )
        // fake initFactoryPtr: no default constructor, only for restart
        .def( py::init( &initFactoryPtr< MeshCoordinatesField, std::string > ) )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( py::self * float() )
        .def( float() * py::self )
        .def( py::self *= float() )
        .def( -py::self )
        .def( "__add__",
              +[]( MeshCoordinatesField &a, const FieldOnNodesReal &b ) { return a + b; } )
        .def( "__add__",
              +[]( const FieldOnNodesReal &a, MeshCoordinatesField &b ) { return a + b; } )
        .def( "__getitem__",
              +[]( const MeshCoordinatesField &v, int i ) { return v.operator[]( i ); }, R"(
Return the coordinate at index *idx* in the vector.

The value is the same as *getValues()[idx]* without creating the entire vector.

Returns:
    float: Values of the *idx*-th coordinate.
        )",
              py::arg( "idx" ) )
        .def( "getValues", &MeshCoordinatesField::getValues, R"(
Return a list of values of the coordinates as (x1, y1, z1, x2, y2, z2...)

Returns:
    list[float]: List of coordinates (size = 3 * number of nodes).
        )" )
        .def( "duplicate", &MeshCoordinatesField::duplicate, R"(
Return a copy of MeshCoordinatesField object

Returns:
    MeshCoordinatesField : MeshCoordinatesField object
        )" )
        .def( "size", &MeshCoordinatesField::size, R"(
Return the size of the field

Returns:
    int : number of values of MeshCoordinatesField object
        )" )
        .def( "updateValuePointers", &MeshCoordinatesField::updateValuePointers, R"(
Update values of internal pointer.
        )" );
};
