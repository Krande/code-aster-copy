/**
 * @file FunctionInterface.cxx
 * @brief Interface python de Function
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

#include "PythonBindings/FunctionInterface.h"

#include "aster_pybind.h"

void exportFunctionToPython( py::module_ &mod ) {

    py::class_< BaseFunction, BaseFunction::BaseFunctionPtr, GenericFunction >( mod,
                                                                                "BaseFunction" )
        // fake initFactoryPtr: created by subclasses
        // fake initFactoryPtr: created by subclasses
        .def( "setParameterName", &Function::setParameterName )
        .def( "setResultName", &Function::setResultName )
        .def( "setInterpolation", &Function::setInterpolation )
        .def( "setValues", &Function::setValues )
        .def( "getValues", &Function::getValues, R"(
Return a list of the values of the function as (x1, x2, ..., y1, y2, ...)

Returns:
    list[float]: List of values (size = 2 * *size()*).

        )" );

    py::class_< Function, Function::FunctionPtr, BaseFunction >( mod, "Function" )
        .def( py::init( &initFactoryPtr< Function > ) )
        .def( py::init( &initFactoryPtr< Function, std::string > ) )
        .def( "setValues", &Function::setValues )
        .def( "size", &Function::size, R"(
Return the number of points of the function.

Returns:
    int: Number of points.

        )" )
        .def( "setAsConstant", &Function::setAsConstant );

    py::class_< FunctionComplex, FunctionComplex::FunctionComplexPtr, BaseFunction >(
        mod, "FunctionComplex" )
        .def( py::init( &initFactoryPtr< FunctionComplex > ) )
        .def( py::init( &initFactoryPtr< FunctionComplex, std::string > ) )
        .def( "setValues", py::overload_cast< const VectorReal &, const VectorReal & >(
                               &FunctionComplex::setValues ) )
        .def( "setValues", py::overload_cast< const VectorReal &, const VectorComplex & >(
                               &FunctionComplex::setValues ) )
        .def( "size", &FunctionComplex::size, R"(
Return the number of points of the function.

Returns:
    int: Number of points.

        )" );
};
