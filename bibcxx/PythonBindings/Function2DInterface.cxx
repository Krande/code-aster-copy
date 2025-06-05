/**
 * @file Function2DInterface.cxx
 * @brief Interface python de Function2D
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

#include "PythonBindings/Function2DInterface.h"

#include "aster_pybind.h"

void exportFunction2DToPython( py::module_ &mod ) {

    py::class_< Function2D, Function2D::Function2DPtr, GenericFunction >( mod, "Function2D" )
        .def( py::init( &initFactoryPtr< Function2D > ) )
        .def( py::init( &initFactoryPtr< Function2D, std::string > ) )
        .def( "getParameters", &Function2D::getParameters, R"(
Return a list of the values of the parameter as (x1, x2, ...)

Returns:
    list[float]: List of values (size = number of functions).

        )" )
        .def( "getValues", &Function2D::getValues, R"(
Return a list of the values of the functions as [F1, F2, ...]
where Fi is (x1, x2, ..., y1, y2, ...).

Returns:
    list[list [float]]: List of values (size = number of functions).
        )" )
        .def( "getProperties", &Function2D::getProperties, R"(
Returns the properties of the function.

Returns:
    tuple[str]: Tuple containing: type of the function (same as `getType()`),
    type of interpolation, parameter name, result name,
    type of extrapolation, object name (same as `getName()`),
    parameter name of functions + a list of dict for each functions that contain
    the type of interpolation and extrapolation.
        )" )
        .def_property( "_t_nappe", &Function2D::getTNappe, &Function2D::setTNappe );
};
