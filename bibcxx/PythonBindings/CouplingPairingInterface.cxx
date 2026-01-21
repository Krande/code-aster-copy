/**
 * @file CouplingPairingInterface.cxx
 * @brief Interface python de CouplingPairing
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
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

#include "PythonBindings/CouplingPairingInterface.h"

#include "aster_pybind.h"
// aslint: disable=C3006

void exportCouplingPairingToPython( py::module_ &mod ) {

    py::enum_< CouplingMethod >( mod, "CouplingMethod", R"(
Enumeration for coupling method.
    )" )
        .value( "Undefined", CouplingMethod::Undefined )
        .value( "Nitsche", CouplingMethod::Nitsche )
        .value( "Penalization", CouplingMethod::Penalization )
        .export_values();

    py::class_< CouplingPairing, CouplingPairingPtr, DataStructure > class_( mod, "CouplingPairing",
                                                                             R"(
Object to create contact pairing.)" );
    class_.def( py::init( &initFactoryPtr< CouplingPairing, ModelPtr, ASTERINTEGER > ) );
    class_.def( "getMesh", &CouplingPairing::getMesh, R"(
Mesh

Returns:
    BaseMesh: the mesh
)" );
    class_.def( "compute", &CouplingPairing::compute,
                R"(
Compute the pairing quantitie

Returns:
    bool: True if the pairing quantities are updated appropriately
)" );
    class_.def( "getNumberOfPairs", &CouplingPairing::getNumberOfPairs, R"(
Return number of pairs

Returns:
    int: number of pairs
)" );
    class_.def( "setVerbosity", &CouplingPairing::setVerbosity, R"(
Set level of verbosity
      0 - without
      1 - normal (default)
      2 - detailled (text)

Arguments:
    level (integer): level of verbosity
        )",
                py::arg( "verbosity" ) );
    class_.def( "getVerbosity", &CouplingPairing::getVerbosity, R"(
Get level of verbosity

Returns:
    integer: level of verbosity
        )" );
    class_.def( "getFiniteElementDescriptor", &CouplingPairing::getFiniteElementDescriptor, R"(
Return Finite Element Descriptor for virtual cells from pairing.

Returns:
    FiniteElementDescriptor: finite element for virtual cells
)" );

    class_.def( "setMasterGroupOfCells", &CouplingPairing::setMasterGroupOfCells, R"(
Set name of the master group.

Arguments:
    group_name [str]: name of the master group.
)",
                py::arg( "group_name" ) );
    class_.def( "setSlaveGroupOfCells", &CouplingPairing::setSlaveGroupOfCells, R"(
Set name of the slave group.

Arguments:
    group_name [str]: name of the slave group.
)",
                py::arg( "group_name" ) );
    class_.def( "getPairingField", &CouplingPairing::getPairingField, R"(
Get intersection points

Returns:
    FieldOnCellsReal: intersection points.
)" );
    class_.def( "setMethod", &CouplingPairing::setMethod, R"(
Set method of coupling.

Argument:
    method [CouplingMethod]: Method to use.
)",
                py::arg( "method" ) );
    class_.def( "setCoefficient", &CouplingPairing::setCoefficient, R"(
Set coefficient of penalization.

Argument:
    coeff [float]: coefficient of penalization.
)",
                py::arg( "coeff" ) );
};
