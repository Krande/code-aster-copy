/**
 * @file DOFNumberingInterface.cxx
 * @brief Interface python de DOFNumbering
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
/* person_in_charge: nicolas.sellenet at edf.fr */

#include "PythonBindings/BaseDOFNumberingInterface.h"

#include "aster_pybind.h"

void exportBaseDOFNumberingToPython( py::module_ &mod ) {

    bool ( BaseDOFNumbering::*f1 )( const ModelPtr model, const ListOfLoadsPtr listOfLoads,
                                    const bool verbose ) = &BaseDOFNumbering::computeNumbering;
    bool ( BaseDOFNumbering::*f2 )( const std::vector< BaseDOFNumbering::MatrElem > matrix,
                                    const bool verbose ) = &BaseDOFNumbering::computeNumbering;

    py::class_< BaseDOFNumbering, BaseDOFNumbering::BaseDOFNumberingPtr, DataStructure >(
        mod, "BaseDOFNumbering" )
        // fake initFactoryPtr: created by subclasses
        // fake initFactoryPtr: created by subclasses
        .def( "computeNumbering", f1, py::arg( "model" ), py::arg( "listOfLoads" ),
              py::arg( "verbose" ) = true )
        .def( "computeNumbering", f2, py::arg( "matrix" ), py::arg( "verbose" ) = true )
        .def( "computeRenumbering", &BaseDOFNumbering::computeRenumbering, py::arg( "model" ),
              py::arg( "listOfLoads" ), py::arg( "defiCont" ), py::arg( "vContElem" ),
              py::arg( "verbose" ) = true )
        .def( "getFiniteElementDescriptors", &BaseDOFNumbering::getFiniteElementDescriptors, R"(
Returns the objects defining the finite elements.

Returns:
    list[FiniteElementDescriptor]: List of finite elements descriptions.
    )" )
        .def( "setFiniteElementDescriptors", &BaseDOFNumbering::setFiniteElementDescriptors, R"(
Returns the object defining the finite elements.

Arguments:
    descr (list[FiniteElementDescriptor]): List of finite elements descriptions.
    )",
              py::arg( "descr" ) )
        .def( "getEquationNumbering", &BaseDOFNumbering::getEquationNumbering, R"(
Returns the global equation numbering object

Returns:
    EquationNumbering: global equation numbering.
        )" )
        .def( "getPhysicalQuantity", &BaseDOFNumbering::getPhysicalQuantity, R"(
Returns the name of the physical quantity that is numbered.

Returns:
    str: physical quantity name.
        )" )
        .def( "isParallel", &BaseDOFNumbering::isParallel, R"(
The numbering is distributed across MPI processes for High Performance Computing.

Returns:
    bool: *True* if used, *False* otherwise.
        )" )
        .def( "getMesh", &BaseDOFNumbering::getMesh, R"(
Return the mesh

Returns:
    MeshPtr: a pointer to the mesh
        )" )
        .def( "getModel", &BaseDOFNumbering::getModel )
        .def( "setModel", &BaseDOFNumbering::setModel )
        .def( "getMorseStorage", &BaseDOFNumbering::getMorseStorage );
};
