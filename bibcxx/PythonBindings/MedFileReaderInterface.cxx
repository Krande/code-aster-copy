/**
 * @file MedFileReaderInterface.cxx
 * @brief Interface python de MedFileReader
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

#include "PythonBindings/MedFileReaderInterface.h"

#include "aster_pybind.h"

// Not DataStructures
// aslint: disable=C3006

#ifdef ASTER_HAVE_MED
void exportMedFileReaderToPython( py::module_ &mod ) {

    py::enum_< MedFileAccessType >( mod, "MedFileAccessType", "EntityType", R"(
Enumeration med access type.
    )" )
        .value( "MedReadOnly", MedReadOnly )
        .value( "MedReadWrite", MedReadWrite )
        .value( "MedCreate", MedCreate )
        .export_values();

    py::class_< MedFileReader, MedFileReader::MedFileReaderPtr >( mod, "MedFileReader" )
        .def( py::init( &initFactoryPtr< MedFileReader > ) )
        .def( "__pickling_disabled__", disable_pickling< MedFileReader >() )

        .def( "close", &MedFileReader::close, R"(Close med file)" )
        .def( "getField",
              py::overload_cast< const std::string & >( &MedFileReader::getField, py::const_ ), R"(
Get field from name

Arguments:
    name (str): field name

Returns:
    MedField: med field of name name
            )",
              py::arg( "name" ) )
        .def( "getField", py::overload_cast< int >( &MedFileReader::getField, py::const_ ), R"(
Get field from iterator

Arguments:
    iterator (int): field iterator

Returns:
    MedField: med field
            )",
              py::arg( "iterator" ) )
        .def( "getFieldNames", &MedFileReader::getFieldNames, R"(
Get all field names

Returns:
    list: list of field names
            )" )
        .def( "getFieldNumber", &MedFileReader::getFieldNumber, R"(
Get field number in field

Returns:
    int: field number
            )" )
        .def( "getMesh", &MedFileReader::getMesh, R"(
Get mesh from iterator

Arguments:
    iterator (int): iterator on mesh

Returns:
    MedMesh: med mesh
            )",
              py::arg( "iterator" ) )
        .def( "getMeshNumber", &MedFileReader::getMeshNumber, R"(
Get mesh number

Returns:
    int: mesh number
            )" )
        .def( "getProfileNumber", &MedFileReader::getProfileNumber, R"(
Get profile number

Returns:
    int: profile number
            )" )
        .def( "open", &MedFileReader::open,
              R"(
Open med file

Arguments:
    path (Path|str): path to med file
    accessType (MedFileAccessType): med access type

Returns:
    int: return code (0 if open is ok)
            )",
              py::arg( "path" ), py::arg( "accessType" ) )
        .def( "openParallel", &MedFileReader::openParallel,
              R"(
Open med file in parallel

Arguments:
    path (Path|str): path to med file
    accessType (MedFileAccessType): med access type

Returns:
    int: return code (0 if open is ok)
            )",
              py::arg( "path" ), py::arg( "accessType" ) );
};

#endif
