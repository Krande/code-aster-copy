/**
 * @file ContactPairingInterface.cxx
 * @brief Interface python de ContactPairing
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

#include "PythonBindings/ContactPairingInterface.h"

#include "aster_pybind.h"
// aslint: disable=C3006

void exportContactPairingToPython( py::module_ &mod ) {

    py::class_< ContactPairing, ContactPairingPtr, DataStructure > class_( mod, "ContactPairing",
                                                                           R"(
Object to create contact pairing.)" );
    class_.def( py::init( &initFactoryPtr< ContactPairing, std::string, ContactNewPtr > ) );
    class_.def( py::init( &initFactoryPtr< ContactPairing, ContactNewPtr > ) );
    class_.def( "getMesh", &ContactPairing::getMesh, R"(
Mesh

Returns:
    BaseMesh: the mesh
)" );
    class_.def( "getCoordinates", &ContactPairing::getCoordinates, R"(
Coordinates of nodes used for pairing (almost always different from the intial mesh).

Returns:
    MeshCoordinatesField: the coordinates field
)" );
    class_.def( "updateCoordinates", &ContactPairing::updateCoordinates, R"(
Update the mesh coordinates given a displacement field

Arguments:
    disp (FieldOnNodes) : field for displacement

)",
                ( py::arg( "disp" ) ) );
    class_.def( "setCoordinates", &ContactPairing::setCoordinates, R"(
Set the mesh coordinates field

Arguments:
    coordinates (MeshCoordinatesField) : coordinates to use for pairing
)",
                ( py::arg( "coordinates" ) ) );
    class_.def( "compute", py::overload_cast<>( &ContactPairing::compute ),
                R"(
Compute the pairing quantities on all zones

Returns:
    bool: True if the pairing quantities are updated appropriately
)" );
    class_.def( "compute", py::overload_cast< ASTERINTEGER & >( &ContactPairing::compute ),
                R"(
Compute the pairing quantities on a zone

Arguments:
    zone_index(int)

Returns:
    bool: True if the pairing quantities are updated appropriately
)",
                ( py::arg( "zone_index" ) ) );
    class_.def( "getNumberOfZones", &ContactPairing::getNumberOfZones, R"(
Return the number of zones

Returns:
    int: number of zones
)" );
    class_.def( "getNumberOfPairs",
                py::overload_cast<>( &ContactPairing::getNumberOfPairs, py::const_ ), R"(
Return number of pairs on all zones

Returns:
    int: number of pairs
)" );
    class_.def(
        "getNumberOfPairs",
        py::overload_cast< const ASTERINTEGER & >( &ContactPairing::getNumberOfPairs, py::const_ ),
        R"(
Return the number of pairs on a zone

Arguments:
    zone_index(int)
Returns:
    int: number of pairs
)",
        ( py::arg( "zone_index" ) ) );
    class_.def(
        "getListOfPairs",
        py::overload_cast< const ASTERINTEGER & >( &ContactPairing::getListOfPairs, py::const_ ),
        R"(
Get list of contact pairs for a contact zone

Arguments:
    zone_index(int)

Returns:
    list[tuple[int, int]]: list of contact pairs
)",
        ( py::arg( "zone_index" ) ) );
    class_.def( "getListOfPairs",
                py::overload_cast<>( &ContactPairing::getListOfPairs, py::const_ ),
                R"(
Get list of contact pairs on all zones

Returns:
    list[tuple[int, int]]: list of contact pairs
)" );
    class_.def( "getIntersectionPoints",
                py::overload_cast< ASTERINTEGER &, CoordinatesSpace >(
                    &ContactPairing::getIntersectionPoints, py::const_ ),
                R"(
Get the intersection points between master and slave cells

Arguments:
    zone_index(int) : index of zone
    CoordinatesSpace (CoordinatesSpace): space to describe coordinates

Returns:
    list[pair]: list of pair of coordinates of intersection points
)",
                ( py::arg( "zone_index" ) ),
                py::arg( "CoordinatesSpace" ) = CoordinatesSpace::Global );
    class_.def( "getNumberOfIntersectionPoints",
                py::overload_cast< ASTERINTEGER & >( &ContactPairing::getNumberOfIntersectionPoints,
                                                     py::const_ ),
                R"(
Get list of the number of intersection points beetween a master and slave cells.

Arguments:
    zone_index(int) : index of zone

Returns:
    list: list of number of intersection points
)",
                ( py::arg( "zone_index" ) ) );
    class_.def( "clearPairing", py::overload_cast<>( &ContactPairing::clearPairing ),
                R"(
Clean pairing for all zones

Returns:
    bool: true if the pairing quantities are cleared
)" );
    class_.def( "clearPairing",
                py::overload_cast< const ASTERINTEGER & >( &ContactPairing::clearPairing ),
                R"(
Clean pairing for a zone

Arguments:
    zone_index(int) : index of zone

Returns:
    bool: true if the pairing quantities are cleared
)",
                ( py::arg( "zone_index" ) ) );
    class_.def( "setVerbosity", &ContactPairing::setVerbosity, R"(
Set level of verbosity
      0 - without
      1 - normal (default)
      2 - detailled (text)

Arguments:
    level (integer): level of verbosity
        )",
                py::arg( "verbosity" ) );
    class_.def( "getVerbosity", &ContactPairing::getVerbosity, R"(
Get level of verbosity

Returns:
    integer: level of verbosity
        )" );
    class_.def( "getFiniteElementDescriptor", &ContactPairing::getFiniteElementDescriptor, R"(
Return Finite Element Descriptor for virtual cells from pairing.

Returns:
    FiniteElementDescriptor: finite element for virtual cells
)" );
};
