/**
 * @file ElementaryVectorInterface.cxx
 * @brief Interface python de ElementaryVector
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
#include "PythonBindings/ElementaryVectorInterface.h"

void exportElementaryVectorToPython() {

    FieldOnNodesRealPtr ( ElementaryVector::*c10 )( const BaseDOFNumberingPtr & ) =
        &ElementaryVector::assembleWithLoadFunctions;
    FieldOnNodesRealPtr ( ElementaryVector::*c11 )( const BaseDOFNumberingPtr &,
                                                    const ASTERDOUBLE& ) =
        &ElementaryVector::assembleWithLoadFunctions;


    void ( ElementaryVector::*c3 )( const MechanicalLoadRealPtr & ) =
        &ElementaryVector::addLoad;

    py::class_< ElementaryVector, ElementaryVector::ElementaryVectorPtr,
            py::bases< DataStructure > >( "ElementaryVector", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< ElementaryVector >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ElementaryVector, std::string >))
        .def( "addLoad", c3 )
        .def( "assemble", &ElementaryVector::assemble )
        .def( "assembleWithLoadFunctions", c10 )
        .def( "assembleWithLoadFunctions", c11 )
        .def( "setType", &ElementaryVector::setType )
        .def( "setListOfLoads", &ElementaryVector::setListOfLoads )
        .def( "build", &ElementaryVector::build );

    py::class_< ElementaryVectorDisplacementReal,
            ElementaryVectorDisplacementRealPtr, py::bases< ElementaryVector > >
            ( "ElementaryVectorDisplacementReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ElementaryVectorDisplacementReal > ) )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ElementaryVectorDisplacementReal,
                                                std::string >));

    py::class_< ElementaryVectorTemperatureReal,
            ElementaryVectorTemperatureRealPtr, py::bases< ElementaryVector > >
            ( "ElementaryVectorTemperatureReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ElementaryVectorTemperatureReal > ) )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ElementaryVectorTemperatureReal,
                                                std::string >));

    py::class_< ElementaryVectorPressureComplex,
            ElementaryVectorPressureComplexPtr, py::bases< ElementaryVector > >
            ( "ElementaryVectorPressureComplex", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ElementaryVectorPressureComplex > ) )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ElementaryVectorPressureComplex,
                                                std::string >));
};
