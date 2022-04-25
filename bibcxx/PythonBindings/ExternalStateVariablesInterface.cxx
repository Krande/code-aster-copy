/**
 * @file ExternalStateVariablesInterface.cxx
 * @brief Main interface for ExternalStateVariable
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

#include "PythonBindings/ExternalStateVariablesInterface.h"

#include "aster_pybind.h"

#include "Results/TransientResult.h"

void exportExternalStateVariablesToPython( py::module_ &mod ) {
    py::class_< EvolutionParameter, EvolutionParameterPtr >( mod, "EvolutionParameter" )
        .def( py::init(
                  &initFactoryPtr< EvolutionParameter, const TransientResultPtr &, std::string > ),
              R"(
            Constructor of object

            Arguments:
                result (TransientResult): transient result to define external state variable
                fieldName (str): field in transient result to define external state variable

            )",
              py::arg( "result" ), py::arg( "fieldName" ) )
        .def( "setTimeFunction",
              py::overload_cast< const FormulaPtr & >( &EvolutionParameter::setTimeFunction ), R"(
            Set function to shift results

            Arguments:
                formula (Formula): formula
            )",
              py::arg( "formula" ) )
        .def( "setTimeFunction",
              py::overload_cast< const FunctionPtr & >( &EvolutionParameter::setTimeFunction ), R"(
            Set function to shift results

            Arguments:
                function (Function): function
            )",
              py::arg( "function" ) )
        .def( "setRightExtension", &EvolutionParameter::setRightExtension, R"(
            Set type of the extension to the right of the function to shift the results

            Arguments:
                typeExtension (str): type of extension ('CONSTANT', 'EXCLU', 'LINEAIRE')
            )",
              py::arg( "typeExtension" ) )
        .def( "setLeftExtension", &EvolutionParameter::setLeftExtension, R"(
            Set type of the extension to the left of the function to shift the results

            Arguments:
                typeExtension (str): type of extension ('CONSTANT', 'EXCLU', 'LINEAIRE')
            )",
              py::arg( "typeExtension" ) );

    py::class_< ExternalStateVariable, ExternalStateVariablePtr >( mod, "ExternalStateVariable" )
        .def( py::init( &initFactoryPtr< ExternalStateVariable, std::string, BaseMeshPtr > ) )
        .def( py::init(
            &initFactoryPtr< ExternalStateVariable, std::string, BaseMeshPtr, std::string > ) )
        .def( py::init( &initFactoryPtr< ExternalStateVariable, externVarEnumInt, BaseMeshPtr > ) )
        .def( py::init(
            &initFactoryPtr< ExternalStateVariable, externVarEnumInt, BaseMeshPtr, std::string > ) )
        .def( "setField", &ExternalStateVariable::setField, R"(
            Define constant value in time for external state variable

            Arguments:
                field (field): field to define value
            )",
              py::arg( "field" ) )
        .def( "getField", &ExternalStateVariable::getField, R"(
            Get the field of values
            )" )
        .def( "getTransientResult", &ExternalStateVariable::getTransientResult, R"(
            Get the transient result
            )" )
        .def( "setReferenceValue", &ExternalStateVariable::setReferenceValue, R"(
            Set reference value for external state variable

            Arguments:
                value (float): reference value
            )",
              py::arg( "value" ) )
        .def( "setEvolutionParameter", &ExternalStateVariable::setEvolutionParameter, R"(
            Define evolution parameters for values of external state variable

            Arguments:
                evolutionParameter (EvolutionParameter): object EvolutionParameter to define
            )",
              py::arg( "evolutionParameter" ) );
};
