/**
 * @file LinearSolverInterface.cxx
 * @brief Interface python de LinearSolver
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include <boost/python.hpp>

namespace py = boost::python;
#include <PythonBindings/factory.h>
#include "PythonBindings/LinearSolverInterface.h"

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( solveWithDirichletBC_overloads,
                                        solveWithDirichletBC, 3, 4 )
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( solve_overloads, solve,
                                        2, 3 )

void exportLinearSolverToPython() {

    py::enum_< LinearSolverEnum >( "BaseLinearSolverName" )
        .value( "MultFront", MultFront )
        .value( "Ldlt", Ldlt )
        .value( "Mumps", Mumps )
        .value( "Petsc", Petsc )
        .value( "Gcpc", Gcpc );

    py::enum_< Renumbering >( "Renumbering" )
        .value( "MD", MD )
        .value( "MDA", MDA )
        .value( "Metis", Metis )
        .value( "RCMK", RCMK )
        .value( "AMD", AMD )
        .value( "AMF", AMF )
        .value( "PORD", PORD )
        .value( "QAMD", QAMD )
        .value( "Scotch", Scotch )
        .value( "Auto", Auto )
        .value( "Parmetis", Parmetis )
        .value( "Ptscotch", Ptscotch )
        .value( "Sans", Sans );

    py::enum_< Preconditioning >( "Preconditioning" )
        .value( "IncompleteLdlt", IncompleteLdlt )
        .value( "SimplePrecisionLdlt", SimplePrecisionLdlt )
        .value( "Jacobi", Jacobi )
        .value( "Sor", Sor )
        .value( "Ml", Ml )
        .value( "Boomer", Boomer )
        .value( "Gamg", Gamg )
        .value( "LagrBloc", LagrBloc )
        .value( "Without", Without );

    py::enum_< IterativeSolverAlgorithm >( "IterativeSolverAlgorithm" )
        .value( "ConjugateGradiant", ConjugateGradiant )
        .value( "ConjugateResidual", ConjugateResidual )
        .value( "GMRes", GMRes )
        .value( "GCR", GCR )
        .value( "FGMRes", FGMRes );

    py::enum_< LagrangeTreatment >( "LagrangeTreatment" )
        .value( "Eliminate", Eliminate )
        .value( "NotEliminate", NotEliminate )
        .value( "LagrangeEliminateReal", LagrangeEliminateReal );

    py::enum_< MemoryManagement >( "MemoryManagement" )
        .value( "InCore", InCore )
        .value( "OutOfCore", OutOfCore )
        .value( "Automatic", Automatic )
        .value( "Evaluation", Evaluation );

    py::enum_< MatrixType >( "MatrixType" )
        .value( "NonSymetric", NonSymetric )
        .value( "Symetric", Symetric )
        .value( "SymetricPositiveDefinite", SymetricPositiveDefinite )
        .value( "Undefined", Undefined );

    py::enum_< MumpsPostTreatment >( "MumpsPostTreatment" )
        .value( "WithoutPostTreatment", WithoutPostTreatment )
        .value( "AutomaticPostTreatment", AutomaticPostTreatment )
        .value( "ForcedPostTreatment", ForcedPostTreatment )
        .value( "MinimalPostTreatment", MinimalPostTreatment );

    py::enum_< MumpsAcceleration >( "MumpsAcceleration" )
        .value( "AutomaticAcceleration", AutomaticAcceleration )
        .value( "FullRank", FullRank )
        .value( "FullRankPlus", FullRankPlus )
        .value( "LowRank", LowRank )
        .value( "LowRankPlus", LowRankPlus );

    py::class_< BaseLinearSolver, BaseLinearSolver::BaseLinearSolverPtr,
                py::bases< DataStructure > >( "BaseLinearSolver", py::no_init )
        .def( "build", &BaseLinearSolver::build )
        .def( "solve", &BaseLinearSolver::solve,
              solve_overloads() )
        .def( "solveWithDirichletBC",
              &BaseLinearSolver::solveWithDirichletBC,
              solveWithDirichletBC_overloads() )
        .def( "disablePreprocessing", &BaseLinearSolver::disablePreprocessing )
        .def( "factorize", &BaseLinearSolver::factorize )
        .def( "setAlgorithm", &BaseLinearSolver::setAlgorithm )
        .def( "setDistributedMatrix", &BaseLinearSolver::setDistributedMatrix )
        .def( "setErrorOnMatrixSingularity",
              &BaseLinearSolver::setErrorOnMatrixSingularity )
        .def( "setFilling", &BaseLinearSolver::setFilling )
        .def( "setFillingLevel", &BaseLinearSolver::setFillingLevel )
        .def( "setLagrangeElimination", &BaseLinearSolver::setLagrangeElimination )
        .def( "setLowRankSize", &BaseLinearSolver::setLowRankSize )
        .def( "setLowRankThreshold", &BaseLinearSolver::setLowRankThreshold )
        .def( "setMatrixFilter", &BaseLinearSolver::setMatrixFilter )
        .def( "setMatrixType", &BaseLinearSolver::setMatrixType )
        .def( "setMaximumNumberOfIteration",
              &BaseLinearSolver::setMaximumNumberOfIteration )
        .def( "setMemoryManagement", &BaseLinearSolver::setMemoryManagement )
        .def( "setPivotingMemory", &BaseLinearSolver::setPivotingMemory )
        .def( "setPrecisionMix", &BaseLinearSolver::setPrecisionMix )
        .def( "setPreconditioning", &BaseLinearSolver::setPreconditioning )
        .def( "setPreconditioningResidual", &BaseLinearSolver::setPreconditioningResidual )
        .def( "setSingularityDetectionThreshold",
              &BaseLinearSolver::setSingularityDetectionThreshold )
        .def( "setSolverResidual", &BaseLinearSolver::setSolverResidual )
        .def( "setPetscOption", &BaseLinearSolver::setPetscOption )
        .def( "setComponentName", &BaseLinearSolver::setComponentName )
        .def( "setComponentPartition", &BaseLinearSolver::setComponentPartition )
        .def( "setUpdatePreconditioningParameter",
              &BaseLinearSolver::setUpdatePreconditioningParameter );

    py::class_< MultFrontSolver, MultFrontSolverPtr,
                py::bases< BaseLinearSolver > >( "MultFrontSolver", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< MultFrontSolver >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< MultFrontSolver, std::string >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< MultFrontSolver, Renumbering >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< MultFrontSolver, std::string, Renumbering >));

    py::class_< LdltSolver, LdltSolverPtr, py::bases< BaseLinearSolver > >(
        "LdltSolver", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< LdltSolver >))
        .def( "__init__", py::make_constructor(&initFactoryPtr< LdltSolver, std::string >))
        .def( "__init__", py::make_constructor(&initFactoryPtr< LdltSolver, Renumbering >))
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< LdltSolver, std::string, Renumbering >));

    py::class_< MumpsSolver, MumpsSolverPtr, py::bases< BaseLinearSolver > >(
        "MumpsSolver", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< MumpsSolver >))
        .def( "__init__", py::make_constructor(&initFactoryPtr< MumpsSolver, std::string >))
        .def( "__init__", py::make_constructor(&initFactoryPtr< MumpsSolver, Renumbering >))
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< MumpsSolver, std::string, Renumbering >))
        .def( "setAcceleration", &BaseLinearSolver::setAcceleration )
        .def( "setPostTreatment", &BaseLinearSolver::setPostTreatment );

    py::class_< PetscSolver, PetscSolverPtr, py::bases< BaseLinearSolver > >(
        "PetscSolver", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< PetscSolver >))
        .def( "__init__", py::make_constructor(&initFactoryPtr< PetscSolver, std::string >))
        .def( "__init__", py::make_constructor(&initFactoryPtr< PetscSolver, Renumbering >))
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< PetscSolver, std::string, Renumbering >));

    py::class_< GcpcSolver, GcpcSolverPtr, py::bases< BaseLinearSolver > >(
        "GcpcSolver", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< GcpcSolver >))
        .def( "__init__", py::make_constructor(&initFactoryPtr< GcpcSolver, std::string >))
        .def( "__init__", py::make_constructor(&initFactoryPtr< GcpcSolver, Renumbering >))
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< GcpcSolver, std::string, Renumbering >));
};
