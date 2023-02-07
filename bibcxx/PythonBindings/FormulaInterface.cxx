/**
 * @file FormulaInterface.cxx
 * @brief Interface python de Formula
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

/* person_in_charge: mathieu.courtois@edf.fr */

#include "PythonBindings/FormulaInterface.h"

#include "aster_pybind.h"

void exportFormulaToPython( py::module_ &mod ) {

    py::class_< Formula, Formula::FormulaPtr, GenericFunction >( mod, "Formula" )
        .def( py::init( &initFactoryPtr< Formula > ) )
        .def( py::init( &initFactoryPtr< Formula, std::string > ) )

        .def( "setVariables", &Formula::setVariables, R"(
Define the variables names.

Arguments:
    varnames (list[str]): List of variables names.
        )",
              py::arg( "varnames" ) )

        .def( "setExpression", &Formula::setExpression, R"(
Define the expression of the formula.

If the expression uses non builtin objects, the evaluation context must be
defined using `:func:setContext`.

Arguments:
    expression (str): Expression of the formula.
        )",
              py::arg( "expression" ) )

        .def( "setComplex", &Formula::setComplex, R"(
Set the type of the formula as complex.
        )" )

        .def( "setContext", &Formula::setContext, R"(
Define the context holding objects required to evaluate the expression.

Arguments:
    context (dict): Context for the evaluation.
        )",
              py::arg( "context" ) )

        .def( "evaluate", &Formula::evaluate, R"(
Evaluate the formula with the given variables values.

Arguments:
    val (list[float]): List of the values of the variables.

Returns:
    float/complex: Value of the formula for these values.
        )",
              py::arg( "*val" ) )

        .def( "getVariables", &Formula::getVariables, R"(
Return the variables names.

Returns:
    list[str]: List of the names of the variables.
        )" )

        .def( "getExpression", &Formula::getExpression, R"(
Return expression of the formula.

Returns:
    str: Expression of the formula.
        )" )

        .def( "getContext", &Formula::getContext, R"(
Return the context used to evaluate the formula.

Returns:
    dict: Context used for evaluation.
        )" )

        .def( "getProperties", &Formula::getProperties, R"(
Return the properties of the formula (for compatibility with function objects).

Returns:
    list[str]: List of 6 strings as function objects contain.
        )" );
};
