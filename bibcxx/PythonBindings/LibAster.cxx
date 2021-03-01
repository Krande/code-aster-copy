/**
 * @file LibAster.cxx
 * @brief Création de LibAster
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

#include <boost/python.hpp>

#include "aster_init.h"
#include "astercxx.h"

#include "PythonBindings/Fortran.h"
#include "Supervis/Exceptions.h"

// Please keep '*Interface.h' files in alphabetical order to ease merging
#include "PythonBindings/AcousticLoadInterface.h"
#include "PythonBindings/AcousticModeResultInterface.h"
#include "PythonBindings/AssemblyMatrixInterface.h"
#include "PythonBindings/BaseDOFNumberingInterface.h"
#include "PythonBindings/BaseMaterialPropertyInterface.h"
#include "PythonBindings/BaseMeshInterface.h"
#include "PythonBindings/BehaviourDefinitionInterface.h"
#include "PythonBindings/BucklingModeResultInterface.h"
#include "PythonBindings/CombinedFourierResultInterface.h"
#include "PythonBindings/ConnectionMeshInterface.h"
#include "PythonBindings/ConstantFieldOnCellsInterface.h"
#include "PythonBindings/ContactInterface.h"
#include "PythonBindings/ConvertersInterface.h"
#include "PythonBindings/CppToFortranGlossaryInterface.h"
#include "PythonBindings/CrackInterface.h"
#include "PythonBindings/CrackShapeInterface.h"
#include "PythonBindings/CrackTipInterface.h"
#include "PythonBindings/CreateEnthalpyInterface.h"
#include "PythonBindings/CyclicSymmetryModeInterface.h"
#include "PythonBindings/DOFNumberingInterface.h"
#include "PythonBindings/DataFieldInterface.h"
#include "PythonBindings/DataStructureInterface.h"
#include "PythonBindings/DebugInterface.h"
#include "PythonBindings/DeleteTemporaryObjectsInterface.h"
#include "PythonBindings/DiscreteProblemInterface.h"
#include "PythonBindings/DynamicMacroElementInterface.h"
#include "PythonBindings/ElasticFourierResultInterface.h"
#include "PythonBindings/ElasticResultInterface.h"
#include "PythonBindings/ElementaryCharacteristicsInterface.h"
#include "PythonBindings/ElementaryMatrixInterface.h"
#include "PythonBindings/ElementaryVectorInterface.h"
#include "PythonBindings/EmpiricalModeResultInterface.h"
#include "PythonBindings/ExternalVariablesConverterInterface.h"
#include "PythonBindings/ExternalVariablesDefinitionInterface.h"
#include "PythonBindings/ExternalVariablesResultInterface.h"
#include "PythonBindings/FiberGeometryInterface.h"
#include "PythonBindings/FieldOnCellsInterface.h"
#include "PythonBindings/FieldOnNodesInterface.h"
#include "PythonBindings/FiniteElementDescriptorInterface.h"
#include "PythonBindings/FluidStructureInteractionInterface.h"
#include "PythonBindings/FluidStructureModalBasisInterface.h"
#include "PythonBindings/FormulaInterface.h"
#include "PythonBindings/FortranInterface.h"
#include "PythonBindings/FullHarmonicAcousticResultInterface.h"
#include "PythonBindings/FullHarmonicResultInterface.h"
#include "PythonBindings/FullResultInterface.h"
#include "PythonBindings/FullTransientResultInterface.h"
#include "PythonBindings/Function2DInterface.h"
#include "PythonBindings/FunctionInterface.h"
#include "PythonBindings/GeneralizedAssemblyMatrixInterface.h"
#include "PythonBindings/GeneralizedAssemblyVectorInterface.h"
#include "PythonBindings/GeneralizedDOFNumberingInterface.h"
#include "PythonBindings/GeneralizedModeResultInterface.h"
#include "PythonBindings/GeneralizedModelInterface.h"
#include "PythonBindings/GeneralizedResultInterface.h"
#include "PythonBindings/GenericFunctionInterface.h"
#include "PythonBindings/GridInterface.h"
#include "PythonBindings/InterspectralMatrixInterface.h"
#include "PythonBindings/DirichletBCInterface.h"
#include "PythonBindings/LinearSolverInterface.h"
#include "PythonBindings/LinearStaticAnalysisInterface.h"
#include "PythonBindings/ListOfFloatsInterface.h"
#include "PythonBindings/ListOfIntegersInterface.h"
#include "PythonBindings/LoadResultInterface.h"
#include "PythonBindings/MPIInfosInterface.h"
#include "PythonBindings/MaterialFieldBuilderInterface.h"
#include "PythonBindings/MaterialFieldInterface.h"
#include "PythonBindings/MaterialInterface.h"
#include "PythonBindings/MaterialPropertyInterface.h"
#include "PythonBindings/MechanicalLoadInterface.h"
#include "PythonBindings/MeshCoordinatesFieldInterface.h"
#include "PythonBindings/MeshEntitiesInterface.h"
#include "PythonBindings/MeshInterface.h"
#include "PythonBindings/MeshesMappingInterface.h"
#include "PythonBindings/ModalBasisInterface.h"
#include "PythonBindings/ModeResultInterface.h"
#include "PythonBindings/ModelInterface.h"
#include "PythonBindings/MultipleElasticResultInterface.h"
#include "PythonBindings/NonLinearResultInterface.h"
#include "PythonBindings/ParallelDOFNumberingInterface.h"
#include "PythonBindings/ParallelFiniteElementDescriptorInterface.h"
#include "PythonBindings/ParallelMechanicalLoadInterface.h"
#include "PythonBindings/ParallelMeshInterface.h"
#include "PythonBindings/PhysicalQuantityInterface.h"
#include "PythonBindings/PhysicsAndModelingsInterface.h"
#include "PythonBindings/PrestressingCableInterface.h"
#include "PythonBindings/ResultInterface.h"
#include "PythonBindings/ResultNamingInterface.h"
#include "PythonBindings/PetscInitializeFinalize.h"
#include "PythonBindings/SimpleFieldOnCellsInterface.h"
#include "PythonBindings/SimpleFieldOnNodesInterface.h"
#include "PythonBindings/SkeletonInterface.h"
#include "PythonBindings/StaticMacroElementInterface.h"
#include "PythonBindings/StructureInterfaceInterface.h"
#include "PythonBindings/StudyDescriptionInterface.h"
#include "PythonBindings/TableContainerInterface.h"
#include "PythonBindings/TableInterface.h"
#include "PythonBindings/ThermalFourierResultInterface.h"
#include "PythonBindings/ThermalLoadInterface.h"
#include "PythonBindings/ThermalResultInterface.h"
#include "PythonBindings/TimeStepperInterface.h"
#include "PythonBindings/TransientResultInterface.h"
#include "PythonBindings/TurbulentSpectrumInterface.h"
#include "PythonBindings/UnitaryThermalLoadInterface.h"
#include "PythonBindings/VariantModalBasisInterface.h"
#include "PythonBindings/VariantStiffnessMatrixInterface.h"
#include "PythonBindings/XfemCrackInterface.h"
// Please keep '*Interface.h' files in alphabetical order to ease merging

namespace py = boost::python;

struct LibAsterInitializer {
    LibAsterInitializer() { initAsterModules(); };

    ~LibAsterInitializer() { jeveux_finalize(); };
};

BOOST_PYTHON_FUNCTION_OVERLOADS( raiseAsterError_overloads, raiseAsterError, 0, 1 )

BOOST_PYTHON_MODULE( libaster ) {
    // hide c++ signatures
    py::docstring_options doc_options( true, true, false );

    boost::shared_ptr< LibAsterInitializer > libGuard( new LibAsterInitializer() );

    py::class_< LibAsterInitializer, boost::shared_ptr< LibAsterInitializer >, boost::noncopyable >(
        "LibAsterInitializer", py::no_init );

    py::scope().attr( "__libguard" ) = libGuard;

    // Definition of exceptions, thrown from 'Exceptions.cxx'/uexcep
    ErrorPy[ASTER_ERROR] = createPyException( "AsterError" );
    py::register_exception_translator< ErrorCpp< ASTER_ERROR > >( &translateError< ASTER_ERROR > );

    ErrorPy[ASTER_CONVERGENCE_ERROR] =
        createPyException( "ConvergenceError", ErrorPy[ASTER_ERROR] );
    py::register_exception_translator< ErrorCpp< ASTER_CONVERGENCE_ERROR > >(
        &translateError< ASTER_CONVERGENCE_ERROR > );

    ErrorPy[ASTER_INTEGRATION_ERROR] =
        createPyException( "IntegrationError", ErrorPy[ASTER_ERROR] );
    py::register_exception_translator< ErrorCpp< ASTER_INTEGRATION_ERROR > >(
        &translateError< ASTER_INTEGRATION_ERROR > );

    ErrorPy[ASTER_SOLVER_ERROR] = createPyException( "SolverError", ErrorPy[ASTER_ERROR] );
    py::register_exception_translator< ErrorCpp< ASTER_SOLVER_ERROR > >(
        &translateError< ASTER_SOLVER_ERROR > );

    ErrorPy[ASTER_CONTACT_ERROR] = createPyException( "ContactError", ErrorPy[ASTER_ERROR] );
    py::register_exception_translator< ErrorCpp< ASTER_CONTACT_ERROR > >(
        &translateError< ASTER_CONTACT_ERROR > );

    ErrorPy[ASTER_TIMELIMIT_ERROR] = createPyException( "TimeLimitError", ErrorPy[ASTER_ERROR] );
    py::register_exception_translator< ErrorCpp< ASTER_TIMELIMIT_ERROR > >(
        &translateError< ASTER_TIMELIMIT_ERROR > );

    py::def( "raiseAsterError", &raiseAsterError, raiseAsterError_overloads() );

    // do not sort (compilation error)
    exportStiffnessMatrixVariantToPython();
    exportModalBasisVariantToPython();
    exportConverters();
    exportDataStructureToPython();
    exportDebugToPython();
    exportMeshEntitiesToPython();
    exportBaseMeshToPython();
    exportMeshToPython();
    exportDiscreteProblemToPython();
    exportBaseDOFNumberingToPython();
    exportDOFNumberingToPython();
    exportElementaryCharacteristicsToPython();
    exportFiniteElementDescriptorToPython();
    exportFiberGeometryToPython();
    exportDataFieldToPython();
    exportFieldOnCellsToPython();
    exportFieldOnNodesToPython();
    exportConstantFieldOnCellsToPython();
    exportSimpleFieldOnCellsToPython();
    exportSimpleFieldOnNodesToPython();
    exportTableToPython();
    exportTableContainerToPython();
    exportTimeStepperToPython();
    exportGeneralizedDOFNumberingToPython();
    exportFluidStructureInteractionToPython();
    exportTurbulentSpectrumToPython();
    exportGenericFunctionToPython();
    exportFunctionToPython();
    exportFormulaToPython();
    exportFortranToPython();
    exportFunction2DToPython();
    exportContactToPython();
    exportAssemblyMatrixToPython();
    exportElementaryMatrixToPython();
    exportElementaryVectorToPython();
    exportGeneralizedAssemblyMatrixToPython();
    exportGeneralizedAssemblyVectorToPython();
    exportInterspectralMatrixToPython();
    exportLinearSolverToPython();
    exportModalBasisToPython();
    exportStructureInterfaceToPython();
    exportAcousticLoadToPython();
    exportDirichletBCToPython();
    exportMechanicalLoadToPython();
    exportPhysicalQuantityToPython();
    exportThermalLoadToPython();
    exportUnitaryThermalLoadToPython();
    exportBehaviourDefinitionToPython();
    exportMaterialToPython();
    exportBaseMaterialPropertyToPython();
    exportMaterialPropertyToPython();
    exportMaterialFieldToPython();
    exportGridToPython();
    exportMeshesMappingToPython();
    exportSkeletonToPython();
    exportDynamicMacroElementToPython();
    exportStaticMacroElementToPython();
    exportCrackShapeToPython();
    exportCrackTipToPython();
    exportCrackToPython();
    exportGeneralizedModelToPython();
    exportModelToPython();
    exportPhysicsAndModelingsToPython();
    exportPrestressingCableToPython();
    exportXfemCrackToPython();
    exportResultToPython();
    exportTransientResultToPython();
    exportLoadResultToPython();
    exportThermalResultToPython();
    exportCombinedFourierResultToPython();
    exportElasticFourierResultToPython();
    exportThermalFourierResultToPython();
    exportMultipleElasticResultToPython();
    exportNonLinearResultToPython();
    exportLinearStaticAnalysisToPython();
    exportStudyDescriptionToPython();
    exportCppToFortranGlossaryToPython();
    exportCyclicSymmetryModeToPython();
    exportFullResultToPython();
    exportModeResultToPython();
    exportModeResultComplexToPython();
    exportAcousticModeResultToPython();
    exportBucklingModeResultToPython();
    exportGeneralizedResultToPython();
    exportElasticResultToPython();
    exportMeshCoordinatesFieldToPython();
    exportFullTransientResultToPython();
    exportFullHarmonicResultToPython();
    exportFullHarmonicAcousticResultToPython();
    exportFluidStructureModalBasisToPython();
    exportGeneralizedModeResultToPython();

#ifdef ASTER_HAVE_MPI
    /* These objects must be declared in Objects/__init__.py as
       OnlyParallelObject for sequential version. */
    exportParallelMeshToPython();
    exportParallelDOFNumberingToPython();
    exportParallelMechanicalLoadToPython();
    exportParallelFiniteElementDescriptorToPython();
#endif /* ASTER_HAVE_MPI */
    exportMPIInfosToPython();

    exportConnectionMeshToPython();
    exportResultNamingToPython();
    exportListOfFloatsToPython();
    exportListOfIntegersToPython();
    exportExternalVariablesToPython();
    exportExternalVariablesConverterToPython();
    exportEmpiricalModeResultToPython();
    exportExternalVariablesResultToPython();
    exportMaterialFieldBuilderToPython();
    exportCreateEnthalpyToPython();
    exportPetscInitializeFinalizeToPython();
    exportDeleteTemporaryObjectsToPython();
};
