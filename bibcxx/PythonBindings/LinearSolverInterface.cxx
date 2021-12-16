/**
 * @file LinearSolverInterface.cxx
 * @brief Interface python de LinearSolver
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include <boost/python.hpp>

namespace py = boost::python;
#include "PythonBindings/LinearSolverInterface.h"
#include <PythonBindings/factory.h>

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( solveWithDirichletBC_overloads, solveWithDirichletBC, 3, 4 )
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( solve_overloads, solve, 2, 3 )

void exportLinearSolverToPython() {

    py::class_< LinearSolver, LinearSolver::LinearSolverPtr, py::bases< DataStructure > >(
        "LinearSolver", py::no_init )
        // .def( "__init__", py::make_constructor( &initFactoryPtr< LinearSolver > ) )
        // .def( "__init__", py::make_constructor( &initFactoryPtr< LinearSolver, std::string >
        // ) )
        .def( "getSolverName", &LinearSolver::getSolverName )
        .def( "supportParallelMesh", &LinearSolver::supportParallelMesh )
        .def( "setKeywords", &LinearSolver::setKeywords )
        .def( "build", &LinearSolver::build )
        .def( "solve", &LinearSolver::solve, solve_overloads() )
        .def( "solveWithDirichletBC", &LinearSolver::solveWithDirichletBC,
              solveWithDirichletBC_overloads() )
        .def( "factorize", &LinearSolver::factorize );

    py::class_< MultFrontSolver, MultFrontSolverPtr, py::bases< LinearSolver > >( "MultFrontSolver",
                                                                                  py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< MultFrontSolver > ) )
        .def( "__init__", py::make_constructor( &initFactoryPtr< MultFrontSolver, std::string > ) );

    py::class_< LdltSolver, LdltSolverPtr, py::bases< LinearSolver > >( "LdltSolver", py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< LdltSolver > ) )
        .def( "__init__", py::make_constructor( &initFactoryPtr< LdltSolver, std::string > ) );

    py::class_< MumpsSolver, MumpsSolverPtr, py::bases< LinearSolver > >( "MumpsSolver",
                                                                          py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< MumpsSolver > ) )
        .def( "__init__", py::make_constructor( &initFactoryPtr< MumpsSolver, std::string > ) );

    py::class_< PetscSolver, PetscSolverPtr, py::bases< LinearSolver > >( "PetscSolver",
                                                                          py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< PetscSolver > ) )
        .def( "__init__", py::make_constructor( &initFactoryPtr< PetscSolver, std::string > ) );

    py::class_< GcpcSolver, GcpcSolverPtr, py::bases< LinearSolver > >( "GcpcSolver", py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< GcpcSolver > ) )
        .def( "__init__", py::make_constructor( &initFactoryPtr< GcpcSolver, std::string > ) );
};
