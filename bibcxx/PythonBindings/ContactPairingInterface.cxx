/**
 * @file ContactPairingInterface.cxx
 * @brief Interface python de ContactPairing
 * @author Nicolas Sellenet
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

#include "PythonBindings/ContactPairingInterface.h"
#include "PythonBindings/factory.h"
#include <boost/python.hpp>

namespace py = boost::python;

void exportContactPairingToPython() {


    py::class_< ContactPairing, ContactPairingPtr, 
                py::bases< DataStructure > >( "ContactPairing", py::no_init )
    .def( "__init__", 
                py::make_constructor( &initFactoryPtr< ContactPairing, 
                            std::vector< ContactZonePtr >, BaseMeshPtr >))
    .def( "__init__", py::make_constructor(
                    &initFactoryPtr< ContactPairing, std::string, 
                            std::vector< ContactZonePtr >, BaseMeshPtr >))
    .def( "getCoordinates", &ContactPairing::getCoordinates, R"(
Compute the new coordinates 
Returns:
    MeshCoordinatesFieldPtr: the new MeshCoordinatesField object
)", ( py::arg( "self" ) ) )
    .def( "updateCoordinates", &ContactPairing::updateCoordinates, R"(
Update the mesh coordinates given a displacement field
)", ( py::arg( "self" ),  py::arg( "disp" ) ) )
    .def( "computePairingQuantities", &ContactPairing::computePairingQuantities, R"(
Compute the pairing quantities associated with the zone izone
Returns:
    bool: True if the pairing quantities are updated appropriately
)", ( py::arg( "self" ), py::arg( "izone" ) ) )
   .def( "getListOfPairs", &ContactPairing::getListOfPairs,  R"(
return list of pairs associated with the zone izone
Returns:
    List[int]: List of pairs
)", ( py::arg( "self" ), py::arg( "izone" ) ) )
   .def( "getNumberOfPairs", &ContactPairing::getNumberOfPairs,  R"(
return number of  pairs associated with the zone izone
Returns:
    int: number of pairs
)", ( py::arg( "self" ), py::arg( "izone" ) ) );

};
