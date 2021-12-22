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

void exportLinearSolverToPython() {

    bool ( LinearSolver::*c1 )( const AssemblyMatrixDisplacementRealPtr ) =
        &LinearSolver::factorize;
    bool ( LinearSolver::*c2 )( const AssemblyMatrixDisplacementComplexPtr ) =
        &LinearSolver::factorize;
    bool ( LinearSolver::*c3 )( const AssemblyMatrixTemperatureRealPtr ) = &LinearSolver::factorize;
    bool ( LinearSolver::*c4 )( const AssemblyMatrixPressureRealPtr ) = &LinearSolver::factorize;

    FieldOnNodesRealPtr ( LinearSolver::*c5 )( const FieldOnNodesRealPtr currentRHS ) const =
        &LinearSolver::solve;
    FieldOnNodesComplexPtr ( LinearSolver::*c6 )( const FieldOnNodesComplexPtr currentRHS ) const =
        &LinearSolver::solve;

    FieldOnNodesRealPtr ( LinearSolver::*c7 )( const FieldOnNodesRealPtr currentRHS,
                                               const FieldOnNodesRealPtr dirichletBCField ) const =
        &LinearSolver::solve;
    FieldOnNodesComplexPtr ( LinearSolver::*c8 )( const FieldOnNodesComplexPtr currentRHS,
                                                  const FieldOnNodesComplexPtr dirichletBCField )
        const = &LinearSolver::solve;

    py::class_< LinearSolver, LinearSolver::LinearSolverPtr, py::bases< DataStructure > >(
        "LinearSolver", py::no_init )
        // .def( "__init__", py::make_constructor( &initFactoryPtr< LinearSolver > ) )
        // .def( "__init__", py::make_constructor( &initFactoryPtr< LinearSolver, std::string >
        // ) )
        .def( "getSolverName", &LinearSolver::getSolverName )
        .def( "supportParallelMesh", &LinearSolver::supportParallelMesh )
        .def( "setKeywords", &LinearSolver::setKeywords )
        .def( "setCommandName", &LinearSolver::setCommandName )
        .def( "enableXfem", &LinearSolver::enableXfem )
        .def( "build", &LinearSolver::build )
        .def( "solve", c5 )
        .def( "solve", c6 )
        .def( "solve", c7 )
        .def( "solve", c8 )
        .def( "factorize", c1 )
        .def( "factorize", c2 )
        .def( "factorize", c3 )
        .def( "factorize", c4 )
        .def( "deleteFactorizedMatrix", &LinearSolver::deleteFactorizedMatrix );

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
