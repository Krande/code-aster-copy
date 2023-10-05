/**
 * @file MedFieldInterface.cxx
 * @brief Interface python de MedField
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

#include "PythonBindings/MedFieldInterface.h"

#include "aster_pybind.h"

// Not DataStructures
// aslint: disable=C3006

#ifdef ASTER_HAVE_MED
void exportMedFieldToPython( py::module_ &mod ) {

    py::class_< MedField, MedField::MedFieldPtr >( mod, "MedField" )
        .def( "__pickling_disabled__", disable_pickling< MedField >() )

        .def( "getComponentName", &MedField::getComponentName, R"(
Get field component name )" )
        .def( "getComponentNumber", &MedField::getComponentNumber, R"(
Get field component number )" )
        .def( "getName", &MedField::getName, R"(
Get field name )" )
        .def( "getAllSupportEntitiesAtSequence", &MedField::getAllSupportEntitiesAtSequence, R"(
Get list of all entity type and geometric type in calculation sequence

Arguments:
    numdt (int): time step id
    numit (int): iteration id

Returns:
    list: list of pair of entity type and geometry type
            )",
              py::arg( "numdt" ), py::arg( "numit" ) )
        .def( "getProfileNumberAtSequenceOnEntity", &MedField::getProfileNumberAtSequenceOnEntity,
              R"(
Get profile number in calculation sequence for a given entity and geometric type )" )
        .def( "getSequence", &MedField::getSequence, R"(
Get time step id and iteration id for a given sequence id

Returns:
    list: time step id and iteration id
            )" )
        .def( "getTime", &MedField::getTime, R"(
Get time for a given sequence id

Returns:
    float: time 
            )" )
        .def( "getSequenceNumber", &MedField::getSequenceNumber, R"(
Get calculation sequence number )" )
        .def( "getValuesAtSequenceOnCellTypesList", &MedField::getValuesAtSequenceOnCellTypesList,
              R"(
Get cell field values at calculation sequence from geometric type list

Arguments:
    numdt (int): time step id
    numit (int): iteration id
    geomtyp (list): list of geomtric types

Returns:
    list: values on cells (same sort as list of geomtric types)
            )",
              py::arg( "numdt" ), py::arg( "numit" ), py::arg( "geomtyp" ) )
        .def( "getValuesAtSequenceOnNodes", &MedField::getValuesAtSequenceOnNodes, R"(
Get node field values at calculation sequence

Arguments:
    numdt (int): time step id
    numit (int): iteration id

Returns:
    list: values on nodes
            )",
              py::arg( "numdt" ), py::arg( "numit" ) )
        .def( "getValuesAtSequenceOnEntityAndProfile",
              &MedField::getValuesAtSequenceOnEntityAndProfile, R"(
Get field values

Arguments:
    numdt (int): time step id
    numit (int): iteration id
    entity (int): entity type
    geometry (int): geometric type
    iterator (int): iterator on profile

Returns:
    list: values
            )",
              py::arg( "numdt" ), py::arg( "numit" ), py::arg( "entity" ), py::arg( "geometry" ),
              py::arg( "iterator" ) );
};

#endif
