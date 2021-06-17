/**
 * @file GeneralizedModeResultInterface.cxx
 * @brief Interface python de GeneralizedModeResult
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

#include "PythonBindings/GeneralizedModeResultInterface.h"
#include "PythonBindings/VariantStiffnessMatrixInterface.h"
#include "PythonBindings/factory.h"
#include <boost/python.hpp>

namespace py = boost::python;

void exportGeneralizedModeResultToPython() {

    bool ( GeneralizedModeResult::*c1 )( const GeneralizedAssemblyMatrixRealPtr & ) =
        &GeneralizedModeResult::setStiffnessMatrix;
    bool ( GeneralizedModeResult::*c2 )( const GeneralizedAssemblyMatrixComplexPtr & ) =
        &GeneralizedModeResult::setStiffnessMatrix;

    py::class_< GeneralizedModeResult, GeneralizedModeResultPtr,
                py::bases< FullResult > >( "GeneralizedModeResult",
                                                             py::no_init )
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< GeneralizedModeResult, std::string >))
        .def( "__init__", py::make_constructor(&initFactoryPtr< GeneralizedModeResult >))
        .def( "setDampingMatrix", &GeneralizedModeResult::setDampingMatrix )
        .def( "getGeneralizedDOFNumbering",
              &GeneralizedModeResult::getGeneralizedDOFNumbering )
        .def( "setGeneralizedDOFNumbering",
              &GeneralizedModeResult::setGeneralizedDOFNumbering )
        .def( "setStiffnessMatrix", c1 )
        .def( "setStiffnessMatrix", c2 )
        .def( "getDampingMatrix", &GeneralizedModeResult::getDampingMatrix )
        .def( "getStiffnessMatrix", &getGeneralizedStiffnessMatrix< GeneralizedModeResultPtr > );
};
