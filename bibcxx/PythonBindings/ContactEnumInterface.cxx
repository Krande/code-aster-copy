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

#include <boost/python.hpp>

namespace py = boost::python;
#include "PythonBindings/ContactEnumInterface.h"

void exportContactEnumToPython() {

    py::enum_< ContactAlgo >( "ContactAlgo" )
        .value( "Lagrangian", ContactAlgo::Lagrangian )
        .value( "Nitsche", ContactAlgo::Nitsche )
        .value( "Penalization", ContactAlgo::Penalization );

    py::enum_< ContactVariant >( "ContactVariant" )
        .value( "Empty", ContactVariant::Empty )
        .value( "Rapide", ContactVariant::Rapide )
        .value( "Robust", ContactVariant::Robust );

    py::enum_< ContactType >( "ContactType" )
        .value( "Unilateral", ContactType::Unilateral )
        .value( "Bilateral", ContactType::Bilateral )
        .value( "Stick", ContactType::Stick );

    py::enum_< FrictionAlgo >( "FrictionAlgo" )
        .value( "Lagrangian", FrictionAlgo::Lagrangian )
        .value( "Nitsche", FrictionAlgo::Nitsche )
        .value( "Penalization", FrictionAlgo::Penalization );

    py::enum_< FrictionType >( "FrictionType" )
        .value( "Without", FrictionType::Without )
        .value( "Tresca", FrictionType::Tresca )
        .value( "Coulomb", FrictionType::Coulomb )
        .value( "Stick", FrictionType::Stick );

    py::enum_< PairingAlgo >( "PairingAlgo" )
        .value( "Mortar", PairingAlgo::Mortar );

    py::enum_< InitialState >( "InitialState" )
        .value( "Interpenetrated", InitialState::Interpenetrated )
        .value( "Yes", InitialState::Yes )
        .value( "No", InitialState::No );
};
