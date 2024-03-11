/**
 * @file CrackInterface.cxx
 * @brief Interface python de Crack
 * @author Nicolas Pignet
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

#include "PythonBindings/CrackInterface.h"

#include "aster_pybind.h"

void exportCrackToPython( py::module_ &mod ) {

    py::class_< Crack, Crack::CrackPtr, DataStructure >( mod, "Crack" )
        .def( py::init( &initFactoryPtr< Crack > ) )
        .def( py::init( &initFactoryPtr< Crack, std::string > ) )
        .def( "getCrackFrontNodes", &Crack::getCrackFrontNodes, R"(
            Return the crack front nodes

            Returns:
                list[str]: the crack nodes
        )" )
        .def( "getCrackFrontBasis", &Crack::getCrackFrontBasis, R"(
            Return the crack front basis

            Returns:
                list[float]: the crack front basis
        )" )
        .def( "getCrackFrontPosition", &Crack::getCrackFrontPosition, R"(
            Return the crack front Position

            Returns:
                list[float]: the crack front Position
        )" )
        .def( "getCrackFrontNodeBasis", &Crack::getCrackFrontNodeBasis )
        .def( "getCrackTipCellsType", &Crack::getCrackTipCellsType )
        .def( "getLowerLipGroupName", &Crack::getLowerLipGroupName )
        .def( "getUpperLipGroupName", &Crack::getUpperLipGroupName )
        .def( "isSymmetric", &Crack::isSymmetric, R"(
            Return true if crack is symeric
        )" )
        .def( "getConfigInit", &Crack::getConfigInit, R"(
            Return the crack initial configuration
        )" );
};
