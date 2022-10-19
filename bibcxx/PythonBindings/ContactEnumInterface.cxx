/**
 * @file ContactEnumInterface.cxx
 * @brief Interface python de ContactEnum
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

#include "PythonBindings/ContactEnumInterface.h"

#include "aster_pybind.h"

void exportContactEnumToPython( py::module_ &mod ) {

    py::enum_< ContactAlgo >( mod, "ContactAlgo" )
        .value( "Lagrangian", ContactAlgo::Lagrangian )
        .value( "Nitsche", ContactAlgo::Nitsche )
        .value( "Penalization", ContactAlgo::Penalization )
        .export_values();

    py::enum_< ContactVariant >( mod, "ContactVariant" )
        .value( "Empty", ContactVariant::Empty )
        .value( "Rapide", ContactVariant::Rapide )
        .value( "Robust", ContactVariant::Robust )
        .value( "Symetric", ContactVariant::Symetric )
        .export_values();

    py::enum_< ContactType >( mod, "ContactType" )
        .value( "Unilateral", ContactType::Unilateral )
        .value( "Bilateral", ContactType::Bilateral )
        .export_values();

    py::enum_< FrictionAlgo >( mod, "FrictionAlgo" )
        .value( "Lagrangian", FrictionAlgo::Lagrangian )
        .value( "Nitsche", FrictionAlgo::Nitsche )
        .value( "Penalization", FrictionAlgo::Penalization )
        .export_values();

    py::enum_< FrictionType >( mod, "FrictionType" )
        .value( "Without", FrictionType::Without )
        .value( "Tresca", FrictionType::Tresca )
        .value( "Coulomb", FrictionType::Coulomb )
        .value( "Stick", FrictionType::Stick )
        .export_values();

    py::enum_< PairingAlgo >( mod, "PairingAlgo" )
        .value( "Mortar", PairingAlgo::Mortar )
        .export_values();

    py::enum_< InitialState >( mod, "InitialState" )
        .value( "Interpenetrated", InitialState::Interpenetrated )
        .value( "Yes", InitialState::Yes )
        .value( "No", InitialState::No )
        .export_values();
};
