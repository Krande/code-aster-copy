/**
 * @file CrackInterface.cxx
 * @brief Interface python de Crack
 * @author Nicolas Pignet
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
        .def( "getCrackFrontAbsCurv", &Crack::getCrackFrontAbsCurv, R"(
            Return the crack front absc curv

            Returns:
                list[float]: the crack front absc curv
        )" )
        .def( "getCrackFrontNodeBasis", &Crack::getCrackFrontNodeBasis, R"(
            Return the basis at each crack front node

            Returns:
                FieldOnNodesReal: field of the crack front basis
        )" )
        .def( "getCrackFrontRadius", &Crack::getCrackFrontRadius, R"(
            Return the crack front Radius

            Returns:
                float: the crack front Radius
        )" )
        .def( "getCrackTipCellsType", &Crack::getCrackTipCellsType, R"(
            Return the crack front cell type

            Returns:
                str: the crack front cell type
        )" )
        .def( "getLowerLipGroupName", &Crack::getLowerLipGroupName, R"(
            Return the group name used to define lower side of cracklip

            Returns:
                str: group name
        )" )
        .def( "getUpperLipGroupName", &Crack::getUpperLipGroupName, R"(
            Return the group name used to define upper side of cracklip

            Returns:
                str: group name
        )" )
        .def( "getLowerNormNodes", &Crack::getLowerNormNodes, R"(
            Return the names for nodes on the lower side of cracklip

            Returns:
                list[str]: node names
        )" )
        .def( "getUpperNormNodes", &Crack::getUpperNormNodes, R"(
            Return the names for nodes on the upper side of cracklip

            Returns:
                list[str]: node names
        )" )
        .def( "getLowerNormNodes2", &Crack::getLowerNormNodes2, R"(
            Return the names for nodes on the lower side of cracklip (for POST_JMOD)

            Returns:
                list[str]: node names
        )" )
        .def( "getUpperNormNodes2", &Crack::getUpperNormNodes2, R"(
            Return the names for nodes on the upper side of cracklip (for POST_JMOD)

            Returns:
                list[str]: node names
        )" )
        .def( "getNormal", &Crack::getNormal, R"(
            Return vector normal of the crack surface

            Returns:
                list[float]: normal to the crack surface
        )" )
        .def( "isSymmetric", &Crack::isSymmetric, R"(
            Return true if crack is symeric
        )" )
        .def( "getConfigInit", &Crack::getConfigInit, R"(
            Return the crack initial configuration
        )" );
};
