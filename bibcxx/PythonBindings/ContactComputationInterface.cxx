/**
 * @file ContactComputationInterface.cxx
 * @brief Interface python de ContactComputation
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

#include "PythonBindings/ContactComputationInterface.h"

#include "aster_pybind.h"

// aslint: disable=C3006

void exportContactComputationToPython( py::module_ &mod ) {

    py::class_< ContactComputation, ContactComputationPtr >( mod, "ContactComputation" )
        .def( py::init( &initFactoryPtr< ContactComputation, ContactNewPtr > ) )
        .def( "geometricGap", &ContactComputation::geometricGap, R"(
Compute geometric gap and indicator using projection. The indicator is equal to 0 for
a node with no projection (gap value is Nan) found else 1.

Arguments:
    coordinates (MeshCoordinatesField): (current) coordinates of mesh

Returns:
    FieldOnNodesReal: gap field.
    FieldOnNodesReal: gap indicator.
        )",
              py::arg( "coordinates" ) )
        .def( "contactData", &ContactComputation::contactData, R"(
Compute contact data (MMCHML)
        )");
};
