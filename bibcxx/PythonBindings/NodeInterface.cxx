/**
 * @brief Interface python de MeshCoordinates
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

#include "astercxx.h"

#include "PythonBindings/NodeInterface.h"

#include "aster_pybind.h"

void exportNodeToPython( py::module_ &mod ) {

    py::class_< Node, NodePtr >( mod, "Node" )
        // .def( py::init( &initFactoryPtr < Node, ASTERINTEGER, std::array< ASTERDOUBLE, 3 > ) )
        // fake initFactoryPtr: not a DataStructure
        .def( define_pickling< Node >() )
        .def(
            "__getitem__",
            +[]( const Node &v, const ASTERINTEGER &i ) { return v.operator[]( i ); } )
        .def(
            "__setitem__",
            +[]( Node &v, const ASTERINTEGER &i, ASTERDOUBLE f ) { return v.operator[]( i ) = f; } )
        .def( "getId", &Node::getId, R"(
Return the Id of the node.

Returns:
    int: local id of the node.
        )" )
        .def( "getValues", &Node::getValues, R"(
Return coordinates as (x,y,z.)

Returns:
    list[float]: (x,y,z).
        )" )
        .def( "x", &Node::x, R"(
Return coordinate x.

Returns:
    float: x.
        )" )
        .def( "y", &Node::y, R"(
Return coordinate y.

Returns:
    float: y.
        )" )
        .def( "z", &Node::z, R"(
Return coordinate z.

Returns:
    float: z.
        )" );
};
