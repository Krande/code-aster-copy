/**
 * @file ContactZoneInterface.cxx
 * @brief Interface python de ContactZone
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

#include "PythonBindings/MeshPairingInterface.h"

#include "aster_pybind.h"
// aslint: disable=C3006

void exportMeshPairingToPython( py::module_ &mod ) {

    py::class_< MeshPairing, MeshPairingPtr, DSWithCppPickling > class_( mod, "MeshPairing", R"(
Object to create a pairing operator between two meshed surfaces.)" );
    py::enum_< CoordinatesSpace >( mod, "CoordinatesSpace", R"(
Type of coordinates: Slave or Global.)" )
        .value( "Slave", CoordinatesSpace::Slave )
        .value( "Global", CoordinatesSpace::Global )
        .export_values();
    class_.def( py::init( &initFactoryPtr< MeshPairing, std::string > ) );
    class_.def( py::init( &initFactoryPtr< MeshPairing > ) );
    class_.def( py::init( &initFactoryPtr< MeshPairing, const py::tuple & > ) );
    class_.def( define_pickling< MeshPairing >() );
    class_.def( "getNumberOfPairs", &MeshPairing::getNumberOfPairs, R"(
Get number of pairs

Returns:
    integer: number of pairs
    )" );
    class_.def( "getListOfPairs", &MeshPairing::getListOfPairs, R"(
Get pairs

Returns:
    list: pairs (slave-master)
        )" );
    class_.def( "getNumberOfIntersectionPoints",
                py::overload_cast< const ASTERINTEGER & >(
                    &MeshPairing::getNumberOfIntersectionPoints, py::const_ ),
                R"(
Get number of intersection points

Arguments:
    indexPair (integer): index of pair

Returns:
    integer: number of intersection points
        )",
                py::arg( "indexPair" ) );
    class_.def(
        "getIntersectionArea",
        py::overload_cast< const ASTERINTEGER & >( &MeshPairing::getIntersectionArea, py::const_ ),
        R"(
Compute intersection of area

Arguments:
    indexPair (integer): index of pair

Returns:
    double: area of intersection
        )",
        py::arg( "indexPair" ) );
    class_.def( "getNumberOfIntersectionPoints",
                py::overload_cast<>( &MeshPairing::getNumberOfIntersectionPoints, py::const_ ),
                R"(
        Get number of intersection points of all pairs

        Returns:
            list: number of intersection points
                )" );
    class_.def( "setPair", &MeshPairing::setPair, R"(
        Set pair of meshed surfaces

        Arguments:
            groupNameSlav (str): slave's name
            groupNameMast (str): master's name
                )",
                py::arg( "groupNameSlav" ), py::arg( "groupNameMast" ) );
    class_.def( "setMesh", &MeshPairing::setMesh, R"(
        Set Mesh

        Arguments:
            mesh (BaseMesh): support mesh
                )",
                py::arg( "mesh" ) );
    class_.def( "compute", &MeshPairing::compute, R"(
Compute pairing

Arguments:
    dist_pairing (real): tolerance from DIST_RATIO (projection outside cell)
    pair_tole (real): tolerance for pairing
        )",
                py::arg( "dist_pairing" ) = -1.0, py::arg( "pair_tole" ) = 1.E-8 );
    class_.def( "setVerbosity", &MeshPairing::setVerbosity, R"(
Set level of verbosity
      0 - without
      1 - normal (default)
      2 - detailled (text)

Arguments:
    level (integer): level of verbosity
        )",
                py::arg( "level" ) );
    class_.def( "getIntersectionPoints", &MeshPairing::getIntersectionPoints, R"(
Get coordinates of intersection points for a given pair

Arguments:
    indexPair (integer): index of pair
    CoordinatesSpace (CoordinatesSpace): space to describe coordinates

Returns:
    list[list]: coordinates in given space
        )",
                py::arg( "indexPair" ), py::arg( "CoordinatesSpace" ) = CoordinatesSpace::Global );
    class_.def( "getIntersectionArea", &MeshPairing::getIntersectionArea, R"(
Get area of intersection for a given pair

Arguments:
    indexPair (integer): index of pair

Returns:
    real: area of intersection
        )",
                py::arg( "indexPair" ) );

    //     class_.def( "getIntersectionPoints", &MeshPairing::getIntersectionPoints, R"(
    // Get coordinates of intersection points for a given pair in global space

    // Arguments:
    //     indexPair (integer): index of pair

    // Returns:
    //     list: intersection points
    //         )",
    //                 py::arg( "indexPair" ) );
    class_.def( "getQuadraturePoints", &MeshPairing::getQuadraturePoints, R"(
Get coordinates of quadrature points for a given pair in global space

Arguments:
    indexPair (integer): index of pair

Returns:
    list: quadrature points
        )",
                py::arg( "indexPair" ) );
    class_.def( "getVerbosity", &MeshPairing::getVerbosity, R"(
Get level of verbosity

Returns:
    integer: level of verbosity
        )" );
    class_.def( "getMesh", &MeshPairing::getMesh, R"(
Return the mesh

Returns:
    Mesh: mesh.
        )" );
    class_.def( "getMasterCellsFromNode", &MeshPairing::getMasterCellsFromNode, R"(
Get the master cells associated with a node number

Arguments:
    int: node number

Returns:
    list: master cells associated
        )",
                py::arg( "node_number" ) );
    class_.def( "getSlaveCellsFromNode", &MeshPairing::getSlaveCellsFromNode, R"(
Get the slave cells associated with a node number

Arguments:
    int: node number

Returns:
    list: slave cells associated
            )",
                py::arg( "node_number" ) );
    class_.def( "getMasterCellNeighbors", &MeshPairing::getMasterCellNeighbors, R"(
Get the master cells in the neighbor of a given master cell number

Arguments:
    int: master cell number

Returns:
    list: master neighbors cells
            )",
                py::arg( "cell_number" ) );
    class_.def( "setExcludedSlaveGroupOfNodes", &MeshPairing::setExcludedSlaveGroupOfNodes,
                R"(
Set excluded groups of nodes on slave side

Arguments:
    str: excluded groups' names
                )",
                py::arg( "groups" ) );
    class_.def( "setExcludedSlaveGroupOfCells", &MeshPairing::setExcludedSlaveGroupOfCells, R"(
Set excluded groups of cells on slave side

Arguments:
    str: excluded groups' names
                )",
                py::arg( "groups" ) );
    class_.def( "getSlaveCellNeighbors", &MeshPairing::getSlaveCellNeighbors, R"(
Get the slave cells in the neighbor of a given slave cell number

Arguments:
    int: slave cell number

Returns:
    list: slave neighbors cells
            )",
                py::arg( "cell_number" ) );
    class_.def( "checkNormals", &MeshPairing::checkNormals, R"(
Check orientation of normals

Arguments:
    ModelPtr: a pointer to the model

Returns:
    nothing
            )",
                py::arg( "model" ) );
    class_.def( "setMethod", &MeshPairing::setMethod, R"(
Set method of pairing

Arguments:
    method (PairingMethod): method ("OLD", "Fast", "Robust)
        )",
                py::arg( "method" ) );
    py::enum_< PairingMethod >( mod, "PairingMethod", R"(
Type of pairing: Fast, BrutForce and Legacy.)" )
        .value( "Fast", PairingMethod::Fast )
        .value( "Legacy", PairingMethod::Legacy )
        .value( "BrutForce", PairingMethod::BrutForce )
        .export_values();
};
