/*
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

#include "aster_pybind.h"

#include "PythonBindings/FieldConverterInterface.h"

void exportFieldConverterToPython( py::module_ &mod ) {

    mod.def( "toFieldOnNodes",
             py::overload_cast< const FieldOnCellsRealPtr >( &toFieldOnNodes< ASTERDOUBLE > ),
             R"(
Convert FieldOnCells to FieldOnNodes

Arguments:
    field (FieldOnCellsReal): field to convert

Returns:
    FieldOnNodesReal: field converted
        )",
             py::arg( "field" ) )
        .def( "toFieldOnNodes",
              py::overload_cast< const FieldOnCellsComplexPtr >( &toFieldOnNodes< ASTERCOMPLEX > ),
              R"(
    Convert FieldOnCells to FieldOnNodes

    Arguments:
        field (FieldOnCellsComplex): field to convert

    Returns:
        FieldOnNodesComplex: field converted
            )",
              py::arg( "field" ) )
        .def(
            "toFieldOnNodes",
            py::overload_cast< const SimpleFieldOnNodesRealPtr >( &toFieldOnNodes< ASTERDOUBLE > ),
            R"(
Convert SimpleFieldOnNodes to FieldOnNodes

Arguments:
    field (SimpleFieldOnNodesReal): field to convert

Returns:
    FieldOnNodesReal: field converted
        )",
            py::arg( "field" ) )
        .def( "toFieldOnNodes",
              py::overload_cast< const SimpleFieldOnNodesComplexPtr >(
                  &toFieldOnNodes< ASTERCOMPLEX > ),
              R"(
    Convert SimpleFieldOnNodes to FieldOnNodes

    Arguments:
        field (SimpleFieldOnNodesComplex): field to convert

    Returns:
        FieldOnNodesComplex: field converted
            )",
              py::arg( "field" ) )
        .def(
            "toSimpleFieldOnNodes",
            py::overload_cast< const FieldOnCellsRealPtr >( &toSimpleFieldOnNodes< ASTERDOUBLE > ),
            R"(
Convert FieldOnCells to SimpleFieldOnNodes

Arguments:
    field (FieldOnCellsReal): field to convert

Returns:
    SimpleFieldOnNodesReal: field converted
        )",
            py::arg( "field" ) )
        .def( "toSimpleFieldOnNodes",
              py::overload_cast< const FieldOnCellsComplexPtr >(
                  &toSimpleFieldOnNodes< ASTERCOMPLEX > ),
              R"(
    Convert FieldOnCells to SimpleFieldOnNodes

    Arguments:
        field (FieldOnCellsComplex): field to convert

    Returns:
        FieldOnNodesComplex: field converted
            )",
              py::arg( "field" ) )
        .def( "toSimpleFieldOnNodes",
              py::overload_cast< const SimpleFieldOnCellsRealPtr >(
                  &toSimpleFieldOnNodes< ASTERDOUBLE > ),
              R"(
Convert SimpleFieldOnCells to SimpleFieldOnNodes

Arguments:
    field (SimpleFieldOnCellsReal): field to convert

Returns:
    SimpleFieldOnNodesReal: field converted
        )",
              py::arg( "field" ) )
        .def(
            "toSimpleFieldOnNodes",
            py::overload_cast< const FieldOnNodesRealPtr >( &toSimpleFieldOnNodes< ASTERDOUBLE > ),
            R"(
Convert FieldOnNodes to SimpleFieldOnNodes

Arguments:
    field (FieldOnNodesReal): field to convert

Returns:
    SimpleFieldOnNodesReal: field converted
        )",
            py::arg( "field" ) )
        .def( "toSimpleFieldOnNodes",
              py::overload_cast< const FieldOnNodesComplexPtr >(
                  &toSimpleFieldOnNodes< ASTERCOMPLEX > ),
              R"(
    Convert FieldOnNodes to SimpleFieldOnNodes

    Arguments:
        field (FieldOnNodesComplex): field to convert

    Returns:
        SimpleFieldOnNodesComplex: field converted
            )",
              py::arg( "field" ) );
};
