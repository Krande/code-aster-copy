/**
 * @file CouplingZonePairingInterface.cxx
 * @brief Interface python de CouplingZonePairing
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

#include "PythonBindings/CouplingZonePairingInterface.h"

#include "aster_pybind.h"
// aslint: disable=C3006

void exportCouplingZonePairingToPython( py::module_ &mod ) {

    py::enum_< CouplingMethod >( mod, "CouplingMethod", R"(
Enumeration for coupling method.
    )" )
        .value( "Undefined", CouplingMethod::Undefined )
        .value( "Nitsche", CouplingMethod::Nitsche )
        .value( "Penalization", CouplingMethod::Penalization )
        .export_values();

    py::class_< CouplingZonePairing, CouplingZonePairingPtr, DataStructure > class_(
        mod, "CouplingZonePairing",
        R"(
Object to create contact pairing.)" );
    class_.def( py::init( &initFactoryPtr< CouplingZonePairing, BaseMeshPtr, ASTERINTEGER > ) );
    class_.def( py::init( &initFactoryPtr< CouplingZonePairing, BaseMeshPtr > ) );

    class_.def( "setMethod", &CouplingZonePairing::setMethod, R"(
Set method.

Returns:
    method [CouplingMethod]: method.
)",
                py::arg( "method" ) );

    class_.def( "setPairingParameters", &CouplingZonePairing::setPairingParameters, R"(
Set pairing parameters.

Arguments:
    parameters [PairingParameter]: PairingParameterPtr.
)",
                py::arg( "parameters" ) );

    class_.def( "setCoefficient", &CouplingZonePairing::setCoefficient, R"(
Set penalization's coefficient.

Arguments:
    coef_pena [float]: penalization's coefficient.
)",
                py::arg( "coef_pena" ) );

    class_.def( "setVerbosity", &CouplingZonePairing::setVerbosity, R"(
Set verbosity.

Arguments:
    verbosity [float]: verbosity level.
)",
                py::arg( "verbosity" ) );

    class_.def( "setSlaveGroupsOfCells", &CouplingZonePairing::setSlaveGroupsOfCells, R"(
Set slave's side.

Arguments:
    groups_name [list[str]]: list of groups.
)",
                py::arg( "groups_name" ) );

    class_.def( "setMasterGroupsOfCells", &CouplingZonePairing::setMasterGroupsOfCells, R"(
Set master's side.

Arguments:
    groups_name [list[str]]: list of groups.
)",
                py::arg( "groups_name" ) );

    class_.def( "check", &CouplingZonePairing::check, R"(
Check common nodes and normals.

Arguments:
    model [Model]: model.
)",
                py::arg( "model" ) );
};
