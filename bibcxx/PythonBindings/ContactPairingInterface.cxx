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

#include "aster_pybind.h"

void exportContactPairingToPython( py::module_ &mod ) {

    py::class_< ContactPairing, ContactPairingPtr, DataStructure >( mod, "ContactPairing" )
        .def( py::init( &initFactoryPtr< ContactPairing, std::string, ContactNewPtr > ) )
        .def( py::init( &initFactoryPtr< ContactPairing, ContactNewPtr > ) )
        .def( "getCoordinates", &ContactPairing::getCoordinates, R"(
Coordinates of nodes used for pairing (almost always different from the intial mesh).

Returns:
    MeshCoordinatesFieldPtr: the coordinates field
)" )
        .def( "updateCoordinates", &ContactPairing::updateCoordinates, R"(
Update the mesh coordinates given a displacement field
)",
              ( py::arg( "disp" ) ) )
        .def( "setCoordinates", &ContactPairing::setCoordinates, R"(
Set the mesh coordinates field

Arguments:
    coordinates (MeshCoordinatesField) : coordinates to use for pairing
)",
              ( py::arg( "coordinates" ) ) )
        .def( "compute", &ContactPairing::compute, R"(
Compute the pairing quantities associated with the zones

Returns:
    bool: True if the pairing quantities are updated appropriately
)" )
        .def( "computeZone", &ContactPairing::computeZone, R"(
Compute the pairing quantities associated with the zone zone_index
Arguments:
    zone_index(int)
Returns:
    bool: True if the pairing quantities are updated appropriately
)",
              ( py::arg( "zone_index" ) ) )
        .def( "getListOfPairsOfZone", &ContactPairing::getListOfPairsOfZone, R"(
return list of pairs associated with the zone izone
Arguments:
    zone_index(int)
Returns:
    List[List[int]]: List of pairs
)",
              ( py::arg( "zone_index" ) ) )
        .def( "getNumberOfPairsOfZone", &ContactPairing::getNumberOfPairsOfZone, R"(
return number of  pairs associated with the zone zone_index
Arguments:
    zone_index(int)
Returns:
    int: number of pairs
)",
              ( py::arg( "zone_index" ) ) )
        .def( "getNumberOfPairs", &ContactPairing::getNumberOfPairs, R"(
return the total number of pairs
Returns:
    int: Total number of pairs
)" )
        .def( "clearZone", &ContactPairing::clearZone, R"(
clean all the paring quantities of zone zonde_index
Arguments:
    zone_index(int)
Returns:
    bool: true if the pairing quantities are cleared
)",
              ( py::arg( "zone_index" ) ) )
        .def( "getSlaveIntersectionPoints",
              py::overload_cast< ASTERINTEGER >( &ContactPairing::getSlaveIntersectionPoints,
                                                 py::const_ ),
              R"(
Get the intersection points beetween a master and slave cells in the parametric
slave space. The maximum number of points is 8.

Arguments:
    zone_index(int) : index of zone

Returns:
    list[list]: list of list of intersection points (each intersection is of size 16)
)",
              ( py::arg( "zone_index" ) ) )
        .def( "clear", &ContactPairing::clear, R"(
clean all the paring quantities of all zones
Returns:
    bool: true if the pairing quantities are cleared
)" )
        .def( "getFiniteElementDescriptor", &ContactPairing::getFiniteElementDescriptor, R"(
Return Finite Element Descriptor for virtual cells from pairing.

Returns:
    FiniteElementDescriptor: finite element for virtual cells
)" );
};
