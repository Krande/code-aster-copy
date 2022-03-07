/**
 * @file ElementaryVectorInterface.cxx
 * @brief Interface python de ElementaryVector
 * @author Nicolas Sellenet
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
#include "PythonBindings/ElementaryVectorInterface.h"
#include "PythonBindings/FieldOnCellsInterface.h"

#include <PythonBindings/factory.h>

void exportElementaryVectorToPython() {

    FieldOnNodesRealPtr ( BaseElementaryVector::*c10 )( const BaseDOFNumberingPtr & ) =
        &BaseElementaryVector::assembleWithLoadFunctions;
    FieldOnNodesRealPtr ( BaseElementaryVector::*c11 )( const BaseDOFNumberingPtr &,
                                                        const ASTERDOUBLE & ) =
        &BaseElementaryVector::assembleWithLoadFunctions;
    FieldOnNodesRealPtr ( BaseElementaryVector::*c12 )( const BaseDOFNumberingPtr &,
                                                        const FieldOnCellsLongPtr &, const int & ) =
        &BaseElementaryVector::assembleWithMask;

    void ( BaseElementaryVector::*c3 )( const MechanicalLoadRealPtr & ) =
        &BaseElementaryVector::addLoad;

    py::class_< BaseElementaryVector, BaseElementaryVectorPtr, py::bases< DataStructure > >(
        "BaseElementaryVector", py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< BaseElementaryVector > ) )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< BaseElementaryVector, std::string > ) )
        .def( "addLoad", c3 )
        .def( "assembleWithLoadFunctions", c10 )
        .def( "assembleWithLoadFunctions", c11 )
        .def( "assembleWithMask", c12 )
        .def( "setType", &BaseElementaryVector::setType )
        .def( "setListOfLoads", &BaseElementaryVector::setListOfLoads )
        .def( "setMaterialField", &BaseElementaryVector::setMaterialField )
        .def( "setModel", &BaseElementaryVector::setModel )
        .def( "setElementaryCharacteristics", &BaseElementaryVector::setElementaryCharacteristics )
        .def( "build", &BaseElementaryVector::build );

    py::class_< ElementaryVectorDisplacementReal, ElementaryVectorDisplacementRealPtr,
                py::bases< BaseElementaryVector > >( "ElementaryVectorDisplacementReal",
                                                     py::no_init )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< ElementaryVectorDisplacementReal > ) )
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< ElementaryVectorDisplacementReal, std::string > ) )
        .def( "assemble", &ElementaryVectorDisplacementReal::assemble );

    py::class_< ElementaryVectorTemperatureReal, ElementaryVectorTemperatureRealPtr,
                py::bases< BaseElementaryVector > >( "ElementaryVectorTemperatureReal",
                                                     py::no_init )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< ElementaryVectorTemperatureReal > ) )
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< ElementaryVectorTemperatureReal, std::string > ) )
        .def( "assemble", &ElementaryVectorTemperatureReal::assemble );

    py::class_< ElementaryVectorPressureComplex, ElementaryVectorPressureComplexPtr,
                py::bases< BaseElementaryVector > >( "ElementaryVectorPressureComplex",
                                                     py::no_init )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< ElementaryVectorPressureComplex > ) )
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< ElementaryVectorPressureComplex, std::string > ) )
        .def( "assemble", &ElementaryVectorPressureComplex::assemble );
};
