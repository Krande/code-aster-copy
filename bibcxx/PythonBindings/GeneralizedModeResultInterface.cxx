/**
 * @file GeneralizedModeResultInterface.cxx
 * @brief Interface python de GeneralizedModeResult
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

#include "PythonBindings/GeneralizedModeResultInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/VariantStiffnessMatrixInterface.h"

void exportGeneralizedModeResultToPython( py::module_ &mod ) {

    py::class_< GeneralizedModeResult, GeneralizedModeResultPtr, FullResult >(
        mod, "GeneralizedModeResult" )
        .def( py::init( &initFactoryPtr< GeneralizedModeResult, std::string > ) )
        .def( py::init( &initFactoryPtr< GeneralizedModeResult > ) )
        .def( "setDampingMatrix", &GeneralizedModeResult::setDampingMatrix )
        .def( "getGeneralizedDOFNumbering", &GeneralizedModeResult::getGeneralizedDOFNumbering )
        .def( "setGeneralizedDOFNumbering", &GeneralizedModeResult::setGeneralizedDOFNumbering )
        .def( "getGeneralizedVectorReal", &GeneralizedModeResult::getGeneralizedVectorReal )
        .def( "getGeneralizedVectorComplex", &GeneralizedModeResult::getGeneralizedVectorComplex )
        .def( "setStiffnessMatrix", py::overload_cast< const GeneralizedAssemblyMatrixRealPtr & >(
                                        &GeneralizedModeResult::setStiffnessMatrix ) )
        .def( "setStiffnessMatrix",
              py::overload_cast< const GeneralizedAssemblyMatrixComplexPtr & >(
                  &GeneralizedModeResult::setStiffnessMatrix ) )
        .def( "getDampingMatrix", &GeneralizedModeResult::getDampingMatrix )
        .def( "getStiffnessMatrix", &getGeneralizedStiffnessMatrix< GeneralizedModeResultPtr > );
};
