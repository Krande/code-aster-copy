/**
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

#include "PythonBindings/PostProcessingInterface.h"

#include "aster_pybind.h"

void exportPostProcessingToPython( py::module_ &mod ) {

    py::class_< PostProcessing, PostProcessing::PostProcessingPtr >( mod, "PostProcessing" )
        .def( py::init( &initFactoryPtr< PostProcessing, PhysicalProblemPtr > ) )
        // fake initFactoryPtr: not a DataStructure
        .def( "computeHydration", &PostProcessing::computeHydration,
              R"(
            Compute hydration at quadrature points (HYDR_ELGA)

            Arguments:
                temp_prev (FieldOnNodesReal): temperature field at begin of current time step
                temp_curr (FieldOnNodesReal): temperature field at end of current time step
                time_prev (float): time at begin of the step
                time_curr (float): time at end of the step
                hydr_prev (FieldOnCellReals): hydration field at begin of current time step

            Returns:
                FieldOnCellReals: hydration field at end of current time step
        )",
              py::arg( "temp_prev" ), py::arg( "temp_curr" ), py::arg( "time_prev" ),
              py::arg( "time_curr" ), py::arg( "hydr_prev" ) )
        .def( "computeAnnealing", &PostProcessing::computeAnnealing,
              R"(
            Modification of internal state variables for annealing

            Arguments:
                internVar (FieldOnNodesReal): internal state variables before annealing
                time_prev (float): time at begin of the step
                time_curr (float): time at end of the step
                externVarPrev (FieldOnCellsReal): external state variables at previous time
                externVarCurr (FieldOnCellsReal): external state variables at current time

            Returns:
                FieldOnCellReals: internal state variables after annealing
        )",
              py::arg( "internVar" ), py::arg( "time_prev" ), py::arg( "time_curr" ),
              py::arg( "externVarPrev" ), py::arg( "externVarCurr" ) )
        .def( "computeMaxResultantForPipe", &PostProcessing::computeMaxResultantForPipe,
              R"(  
            Computes the maximum of the EFGE_ELNO or EGRU_ELNO field in absolute value, 
            based on the maximal values of the equivalent moment at each element.

            Arguments:
             result (Result) : ResultPtr
                The result object containing the fields
             field_name (str) : It should be 'EFGE_ELNO' or 'EGRU_ELNO'

            Returns:
             FieldOnCellReals: The maximal value of the field )",
              py::arg( "result" ), py::arg( "field_name" ) );
};
