/**
 * @file DOFNumberingInterface.cxx
 * @brief Interface python de DOFNumbering
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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
/* person_in_charge: nicolas.sellenet at edf.fr */

#include "PythonBindings/BaseDOFNumberingInterface.h"

#include "aster_pybind.h"

void exportBaseDOFNumberingToPython( py::module_ &mod ) {

    bool ( BaseDOFNumbering::*f1 )( const ModelPtr model, const ListOfLoadsPtr listOfLoads ) =
        &BaseDOFNumbering::computeNumbering;
    bool ( BaseDOFNumbering::*f2 )( const std::vector< BaseDOFNumbering::MatrElem > matrix ) =
        &BaseDOFNumbering::computeNumbering;

    py::class_< BaseDOFNumbering, BaseDOFNumbering::BaseDOFNumberingPtr, DataStructure > c1(
        mod, "BaseDOFNumbering" );
    // fake initFactoryPtr: created by subclasses
    // fake initFactoryPtr: created by subclasses
    c1.def( "addFiniteElementDescriptor", &BaseDOFNumbering::addFiniteElementDescriptor );
    c1.def( "computeNumbering", f1 );
    c1.def( "computeNumbering", f2 );
    c1.def( "computeRenumbering", &BaseDOFNumbering::computeRenumbering );
    c1.def( "getFiniteElementDescriptors", &BaseDOFNumbering::getFiniteElementDescriptors );
    c1.def( "getPhysicalQuantity", &BaseDOFNumbering::getPhysicalQuantity, R"(
Returns the name of the physical quantity that is numbered.

Returns:
    str: physical quantity name.
        )" );

    c1.def( "isParallel", &BaseDOFNumbering::isParallel, R"(
The numbering is distributed across MPI processes for High Performance Computing.

Returns:
    bool: *True* if used, *False* otherwise.
        )" );
    c1.def( "getMesh", &BaseDOFNumbering::getMesh, R"(
Return the mesh

Returns:
    MeshPtr: a pointer to the mesh
        )" );
    c1.def( "getModel", &BaseDOFNumbering::getModel );
    c1.def( "setModel", &BaseDOFNumbering::setModel );
};
