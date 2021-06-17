/**
 * @file TableContainerInterface.cxx
 * @brief Interface python de TableContainer
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

#include <boost/python.hpp>

namespace py = boost::python;
#include <PythonBindings/factory.h>
#include "PythonBindings/TableContainerInterface.h"
#include "PythonBindings/DataStructureInterface.h"

void exportTableContainerToPython() {

    void ( TableContainer::*c1 )( const std::string &,
                                          ElementaryMatrixDisplacementRealPtr ) =
        &TableContainer::addObject;
    void ( TableContainer::*c2 )( const std::string &,
                                          ElementaryMatrixTemperatureRealPtr ) =
        &TableContainer::addObject;
    void ( TableContainer::*c3 )( const std::string &,
                                          ElementaryVectorDisplacementRealPtr ) =
        &TableContainer::addObject;
    void ( TableContainer::*c4 )( const std::string &,
                                          ElementaryVectorTemperatureRealPtr ) =
        &TableContainer::addObject;
    void ( TableContainer::*c5 )( const std::string &, FieldOnCellsRealPtr ) =
        &TableContainer::addObject;
    void ( TableContainer::*c6 )( const std::string &, FieldOnNodesRealPtr ) =
        &TableContainer::addObject;
    void ( TableContainer::*c7 )( const std::string &, FunctionPtr ) =
        &TableContainer::addObject;
    void ( TableContainer::*c8 )( const std::string &, FunctionComplexPtr ) =
        &TableContainer::addObject;
    void ( TableContainer::*c9 )( const std::string &,
                                          GeneralizedAssemblyMatrixRealPtr ) =
        &TableContainer::addObject;
    void ( TableContainer::*c10 )( const std::string &, DataFieldPtr ) =
        &TableContainer::addObject;
    void ( TableContainer::*c11 )( const std::string &, ModeResultPtr ) =
        &TableContainer::addObject;
    void ( TableContainer::*c12 )( const std::string &, ConstantFieldOnCellsRealPtr ) =
        &TableContainer::addObject;
    void ( TableContainer::*c13 )( const std::string &, Function2DPtr ) =
        &TableContainer::addObject;
    void ( TableContainer::*c14 )( const std::string &, TablePtr ) =
        &TableContainer::addObject;

    py::class_< TableContainer, TableContainer::TableContainerPtr,
                py::bases< Table > >( "TableContainer", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< TableContainer >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< TableContainer, std::string >))
        .def( "addObject", c1 )
        .def( "addObject", c2 )
        .def( "addObject", c3 )
        .def( "addObject", c4 )
        .def( "addObject", c5 )
        .def( "addObject", c6 )
        .def( "addObject", c7 )
        .def( "addObject", c8 )
        .def( "addObject", c9 )
        .def( "addObject", c10 )
        .def( "addObject", c11 )
        .def( "addObject", c12 )
        .def( "addObject", c13 )
        .def( "addObject", c14 )
        .def( "getElementaryMatrixDisplacementReal",
              &TableContainer::getElementaryMatrixDisplacementReal )
        .def( "getElementaryMatrixTemperatureReal",
              &TableContainer::getElementaryMatrixTemperatureReal )
        .def( "getElementaryVectorDisplacementReal",
              &TableContainer::getElementaryVectorDisplacementReal )
        .def( "getElementaryVectorTemperatureReal",
              &TableContainer::getElementaryVectorTemperatureReal )
        .def( "getFieldOnCellsReal", &TableContainer::getFieldOnCellsReal )
        .def( "getFieldOnNodesReal", &TableContainer::getFieldOnNodesReal )
        .def( "getFunction", &TableContainer::getFunction )
        .def( "getFunctionComplex", &TableContainer::getFunctionComplex )
        .def( "getGeneralizedAssemblyMatrix",
              &TableContainer::getGeneralizedAssemblyMatrix )
        .def( "getDataField", &TableContainer::getDataField )
        .def( "getModeResult", &TableContainer::getModeResult )
        .def( "getConstantFieldOnCellsReal", &TableContainer::getConstantFieldOnCellsReal )
        .def( "getFunction2D", &TableContainer::getFunction2D )
        .def( "getTable", &TableContainer::getTable );
};
