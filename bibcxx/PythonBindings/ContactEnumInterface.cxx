/**
 * @file ContactEnumInterface.cxx
 * @brief Interface python de ContactEnum
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

#include "PythonBindings/ContactEnumInterface.h"

#include "aster_pybind.h"

void exportContactEnumToPython( py::module_ &mod ) {

    py::enum_< ContactAlgo >( mod, "ContactAlgo", R"(
Enumeration for contact algorithm.
    )" )
        .value( "Lagrangian", ContactAlgo::Lagrangian )
        .value( "Nitsche", ContactAlgo::Nitsche )
        .value( "Penalization", ContactAlgo::Penalization )
        .export_values();

    py::enum_< ContactVariant >( mod, "ContactVariant", R"(
Enumeration for contact variant.
    )" )
        .value( "Empty", ContactVariant::Empty )
        .value( "Fast", ContactVariant::Fast )
        .value( "Robust", ContactVariant::Robust )
        .value( "Symetric", ContactVariant::Symetric )
        .value( "Classic", ContactVariant::Classic )
        .export_values();

    py::enum_< ContactType >( mod, "ContactType", R"(
Enumeration for contact type.
    )" )
        .value( "Unilateral", ContactType::Unilateral )
        .value( "Bilateral", ContactType::Bilateral )
        .export_values();

    py::enum_< FrictionAlgo >( mod, "FrictionAlgo", R"(
Enumeration for friction algorithm.
    )" )
        .value( "Lagrangian", FrictionAlgo::Lagrangian )
        .value( "Nitsche", FrictionAlgo::Nitsche )
        .value( "Penalization", FrictionAlgo::Penalization )
        .export_values();

    py::enum_< FrictionType >( mod, "FrictionType", R"(
Enumeration for friction type.
    )" )
        .value( "Without", FrictionType::Without )
        .value( "Tresca", FrictionType::Tresca )
        .value( "Coulomb", FrictionType::Coulomb )
        .value( "Stick", FrictionType::Stick )
        .export_values();

    py::enum_< PairingAlgo >( mod, "PairingAlgo", R"(
Enumeration for pairing algorithm.
    )" )
        .value( "Mortar", PairingAlgo::Mortar )
        .export_values();

    py::enum_< InitialState >( mod, "InitialState", R"(
Enumeration for initial state.
    )" )
        .value( "Interpenetrated", InitialState::Interpenetrated )
        .value( "Yes", InitialState::Yes )
        .value( "No", InitialState::No )
        .export_values();

    py::enum_< JacobianType >( mod, "JacobianType", R"(
Enumeration for jacobian type.
    )" )
        .value( "Analytical", JacobianType::Analytical )
        .value( "Perturbation", JacobianType::Perturbation )
        .export_values();
};
