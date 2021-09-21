/**
 * @file LinearStaticAnalysisInterface.cxx
 * @brief Interface python de LinearStaticAnalysis
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

#include "PythonBindings/LinearStaticAnalysisInterface.h"
#include "PythonBindings/LoadUtilities.h"
#include "PythonBindings/factory.h"
#include <boost/python.hpp>

namespace py = boost::python;

void exportLinearStaticAnalysisToPython() {

    py::class_< LinearStaticAnalysis, LinearStaticAnalysisPtr > c1(
        "LinearStaticAnalysis", py::no_init );
    c1.def( "__init__",
            py::make_constructor(
                &initFactoryPtr< LinearStaticAnalysis, ModelPtr, MaterialFieldPtr >));
    c1.def( "__init__",
            py::make_constructor(&initFactoryPtr< LinearStaticAnalysis, ModelPtr,
                                              MaterialFieldPtr, ElementaryCharacteristicsPtr >));
    addDirichletBCToInterface( c1 );
    addMechanicalLoadToInterface( c1 );
#ifdef ASTER_HAVE_MPI
    addParallelMechanicalLoadToInterface( c1 );
#endif /* ASTER_HAVE_MPI */
    c1.def( "execute", &LinearStaticAnalysis::execute );
    c1.def( "setLinearSolver", &LinearStaticAnalysis::setLinearSolver );
    c1.def( "setTimeStepManager", &LinearStaticAnalysis::setTimeStepManager );
    c1.def( "setStressComputation", &LinearStaticAnalysis::setStressComputation );
};
