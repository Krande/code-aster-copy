/**
 * @file FunctionInterface.cxx
 * @brief Interface python de Function
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

#include "PythonBindings/FunctionInterface.h"
#include "PythonBindings/factory.h"

void exportFunctionToPython() {

    py::class_< BaseFunction, BaseFunction::BaseFunctionPtr,
                py::bases< GenericFunction > >( "BaseFunction", py::no_init )
        // fake initFactoryPtr: created by subclasses
        // fake initFactoryPtr: created by subclasses
        .def( "setParameterName", &Function::setParameterName )
        .def( "setResultName", &Function::setResultName )
        .def( "setInterpolation", &Function::setInterpolation )
        .def( "setValues", &Function::setValues )
        .def( "getValues", &Function::getValues, R"(
Return a list of the values of the functions as (x1, x2, ..., y1, y2, ...)

Returns:
    list[float]: List of values (size = 2 * *size()*).

        )",
        ( py::arg( "self" ) ) );

    py::class_< Function, Function::FunctionPtr,
                py::bases< BaseFunction > >( "Function", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< Function >))
        .def( "__init__", py::make_constructor(&initFactoryPtr< Function, std::string >))
        .def( "setValues", &Function::setValues )
        .def( "size", &Function::size, R"(
Return the number of points of the function.

Returns:
    int: Number of points.

        )",
        ( py::arg( "self" ) ) )
        .def( "setAsConstant", &Function::setAsConstant );

    // Candidates for setValues
    void ( FunctionComplex::*c1 )( const VectorReal &absc, const VectorReal &ord ) =
        &FunctionComplex::setValues;
    void ( FunctionComplex::*c2 )( const VectorReal &absc, const VectorComplex &ord ) =
        &FunctionComplex::setValues;

    py::class_< FunctionComplex, FunctionComplex::FunctionComplexPtr,
                py::bases< BaseFunction > >( "FunctionComplex", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< FunctionComplex >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< FunctionComplex, std::string >))
        .def( "setValues", c1 )
        .def( "setValues", c2 )
        .def( "size", &FunctionComplex::size, R"(
Return the number of points of the function.

Returns:
    int: Number of points.

        )",
        ( py::arg( "self" ) ) );
};
