/**
 * @file LinearSolverInterface.cxx
 * @brief Interface python de LinearSolver
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

#include "PythonBindings/LinearSolverInterface.h"

#include "aster_pybind.h"

void exportLinearSolverToPython( py::module_ &mod ) {

    py::class_< LinearSolver, LinearSolver::LinearSolverPtr, DataStructure >( mod, "LinearSolver" )
        // .def( py::init( &initFactoryPtr< LinearSolver > ) )
        // .def( py::init( &initFactoryPtr< LinearSolver, std::string >
        // ) )
        .def( "getSolverName", &LinearSolver::getSolverName, R"(
Get the name of the solver used between 'MUMPS', 'PETSC', 'MULT_FRONT' and 'PETSC'

Returns:
     str: name of the solver used
        )" )
        .def( "supportParallelMesh", &LinearSolver::supportParallelMesh, R"(
tell if the solver is enable in HPC

Returns:
     bool: True if the solver support ParallelMesh, else False
        )" )
        .def( "setKeywords", &LinearSolver::setKeywords )
        .def( "setCataPath", &LinearSolver::setCataPath, R"(
Set the path of the catalog that defines the solver keywords.
It can be command name or a path as *code_aster.Cata.Commons.xxxx*.

Arguments:
     path (str): command name or path of the catalog.
        )",
              py::arg( "path" ) )
        .def( "enableXfem", &LinearSolver::enableXfem, R"(
Enable preconditionning for XFEM modeling.
        )" )
        .def( "build", &LinearSolver::build, R"(
build internal objects of the solver

Returns:
     bool: True if the building is a success, else False
        )" )
        .def( "solve",
              py::overload_cast< const FieldOnNodesRealPtr, const FieldOnNodesRealPtr >(
                  &LinearSolver::solve, py::const_ ),
              py::arg( "rhs" ), py::arg( "dirichletBC" ) = nullptr )
        .def( "solve",
              py::overload_cast< const FieldOnNodesComplexPtr, const FieldOnNodesComplexPtr >(
                  &LinearSolver::solve, py::const_ ),
              py::arg( "rhs" ), py::arg( "dirichletBC" ) = nullptr )
        .def( "factorize", &LinearSolver::factorize, R"(
Factorize the matrix.

Arguments:
    matrix (BaseAssemblyMatrix) : matrix to factorize
    raiseException (bool): if *True* an exception is raised in case of error,
    otherwise it stops with an error (default: *False*).

Returns:
    bool: *True* if factorization is a success, else *False*
        )",
              py::arg( "matrix" ), py::arg( "raiseException" ) = false )
        .def( "getMatrix", &LinearSolver::getMatrix, R"(
return the factorized matrix

Returns:
    BaseAssemblyMatrix: factorized matrix
        )" )
        .def( "getPrecondMatrix", &LinearSolver::getPrecondMatrix, R"(
return the preconditionning matrix

Returns:
    BaseAssemblyMatrix: preconditionning matrix
        )" )
        .def( "deleteFactorizedMatrix", &LinearSolver::deleteFactorizedMatrix, R"(
delete the factorized matrix and its preconditionner if created.
This is the case for Mumps and Petsc.

Returns:
     bool: True if success, else False
        )" );

    py::class_< MultFrontSolver, MultFrontSolverPtr, LinearSolver >( mod, "MultFrontSolver" )
        .def( py::init( &initFactoryPtr< MultFrontSolver > ) )
        .def( py::init( &initFactoryPtr< MultFrontSolver, std::string > ) );

    py::class_< LdltSolver, LdltSolverPtr, LinearSolver >( mod, "LdltSolver" )
        .def( py::init( &initFactoryPtr< LdltSolver > ) )
        .def( py::init( &initFactoryPtr< LdltSolver, std::string > ) );

    py::class_< MumpsSolver, MumpsSolverPtr, LinearSolver >( mod, "MumpsSolver" )
        .def( py::init( &initFactoryPtr< MumpsSolver > ) )
        .def( py::init( &initFactoryPtr< MumpsSolver, std::string > ) );

    py::class_< PetscSolver, PetscSolverPtr, LinearSolver >( mod, "PetscSolver" )
        .def( py::init( &initFactoryPtr< PetscSolver > ) )
        .def( py::init( &initFactoryPtr< PetscSolver, std::string > ) )
        .def( "getPetscOptions", &PetscSolver::getPetscOptions, R"(
return the petsc solver options

Returns:
    string: the petsc solver options
        )" );

    py::class_< GcpcSolver, GcpcSolverPtr, LinearSolver >( mod, "GcpcSolver" )
        .def( py::init( &initFactoryPtr< GcpcSolver > ) )
        .def( py::init( &initFactoryPtr< GcpcSolver, std::string > ) );
};
