/**
 * @file BucklingModeResultInterface.cxx
 * @brief Interface python de BucklingModeResult
 * @author Natacha Béreux
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

#include "PythonBindings/BucklingModeResultInterface.h"
#include "PythonBindings/factory.h"
#include <boost/python.hpp>

namespace py = boost::python;
#include "PythonBindings/VariantStiffnessMatrixInterface.h"

void exportBucklingModeResultToPython() {

    bool ( BucklingModeResult::*c1 )( const AssemblyMatrixDisplacementRealPtr & ) =
        &BucklingModeResult::setStiffnessMatrix;
    bool ( BucklingModeResult::*c2 )( const AssemblyMatrixTemperatureRealPtr & ) =
        &BucklingModeResult::setStiffnessMatrix;
    bool ( BucklingModeResult::*c3 )( const AssemblyMatrixDisplacementComplexPtr & ) =
        &BucklingModeResult::setStiffnessMatrix;
    bool ( BucklingModeResult::*c4 )( const AssemblyMatrixPressureRealPtr & ) =
        &BucklingModeResult::setStiffnessMatrix;

    py::class_< BucklingModeResult, BucklingModeResultPtr,
            py::bases< FullResult > >( "BucklingModeResult", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< BucklingModeResult >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< BucklingModeResult, std::string >))
        .def( "getStiffnessMatrix", &getStiffnessMatrix< BucklingModeResultPtr > )
        .def( "setStiffnessMatrix", c1 )
        .def( "setStiffnessMatrix", c2 )
        .def( "setStiffnessMatrix", c3 )
        .def( "setStiffnessMatrix", c4 );
};
