/**
 * @file GenericFunctionInterface.cxx
 * @brief Interface python de GenericFunction
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

#include "PythonBindings/GenericFunctionInterface.h"

#include "aster_pybind.h"

void exportGenericFunctionToPython( py::module_ &mod ) {

    py::class_< GenericFunction, GenericFunction::GenericFunctionPtr, DataStructure >(
        mod, "GenericFunction" )
        // fake initFactoryPtr: created by subclasses
        // fake initFactoryPtr: created by subclasses
        .def( "getProperties", &GenericFunction::getProperties, R"(
Returns the properties of the function.

Returns:
    tuple[str]: Tuple containing: type of the function (same as `getType()`),
    type of interpolation, parameter name, result name,
    type of extrapolation, object name (same as `getName()`).
        )" )
        .def( "setExtrapolation", &GenericFunction::setExtrapolation, R"(
Define the type of extrapolation.

Supported extrapolation types are: "L" for linear, "C" for constant and
"E" for no extrapolation allowed.

Arguments:
    type (str): Type of extrapolation on left and on right. Examples: "CC",
        "LE"...
        )",
              py::arg( "type" ) );
};
