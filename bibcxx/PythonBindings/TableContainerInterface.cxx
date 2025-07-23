/**
 * @file TableContainerInterface.cxx
 * @brief Interface python de TableContainer
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

#include "PythonBindings/TableContainerInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/DataStructureInterface.h"

void exportTableContainerToPython( py::module_ &mod ) {

    py::class_< TableContainer, TableContainer::TableContainerPtr, Table >( mod, "TableContainer" )
        .def( py::init( &initFactoryPtr< TableContainer > ) )
        .def( py::init( &initFactoryPtr< TableContainer, std::string > ) )
        .def( "addObject",
              py::overload_cast< const std::string &, ElementaryMatrixDisplacementRealPtr >(
                  &TableContainer::addObject ) )
        .def( "addObject",
              py::overload_cast< const std::string &, ElementaryVectorDisplacementRealPtr >(
                  &TableContainer::addObject ) )
        .def( "addObject", py::overload_cast< const std::string &, FieldOnCellsRealPtr >(
                               &TableContainer::addObject ) )
        .def( "addObject", py::overload_cast< const std::string &, FieldOnNodesRealPtr >(
                               &TableContainer::addObject ) )
        .def( "addObject", py::overload_cast< const std::string &, FunctionComplexPtr >(
                               &TableContainer::addObject ) )
        .def( "addObject",
              py::overload_cast< const std::string &, GeneralizedAssemblyMatrixRealPtr >(
                  &TableContainer::addObject ) )
        .def( "addObject",
              py::overload_cast< const std::string &, DataFieldPtr >( &TableContainer::addObject ) )
        .def( "addObject", py::overload_cast< const std::string &, ModeResultPtr >(
                               &TableContainer::addObject ) )
        .def( "addObject", py::overload_cast< const std::string &, ConstantFieldOnCellsRealPtr >(
                               &TableContainer::addObject ) )
        .def( "addObject", py::overload_cast< const std::string &, Function2DPtr >(
                               &TableContainer::addObject ) )
        .def( "addObject",
              py::overload_cast< const std::string &, TablePtr >( &TableContainer::addObject ) )
        .def( "addObject",
              py::overload_cast< const std::string &, FunctionPtr >( &TableContainer::addObject ),
              R"(
            Store a *DataStructure* in the table.

            Arguments:
                name (str): String to identify the object in the table.
                object (misc): Object to be inserted, can be a Function, Mesh, Fields...
            )",
              py::arg( "name" ), py::arg( "object" ) )
        .def( "getElementaryMatrixDisplacementReal",
              &TableContainer::getElementaryMatrixDisplacementReal )
        .def( "getElementaryVectorDisplacementReal",
              &TableContainer::getElementaryVectorDisplacementReal )
        .def( "getFieldOnCellsReal", &TableContainer::getFieldOnCellsReal )
        .def( "getFieldOnNodesReal", &TableContainer::getFieldOnNodesReal )
        .def( "getFunction", &TableContainer::getFunction )
        .def( "getFunctionComplex", &TableContainer::getFunctionComplex )
        .def( "getGeneralizedAssemblyMatrix", &TableContainer::getGeneralizedAssemblyMatrix )
        .def( "getDataField", &TableContainer::getDataField )
        .def( "getModeResult", &TableContainer::getModeResult )
        .def( "getConstantFieldOnCellsReal", &TableContainer::getConstantFieldOnCellsReal )
        .def( "getFunction2D", &TableContainer::getFunction2D )
        .def( "getTable", &TableContainer::getTable );
};
