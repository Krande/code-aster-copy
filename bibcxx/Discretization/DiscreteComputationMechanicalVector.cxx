/**
 * @file DiscreteComputation.cxx
 * @brief Implementation of class DiscreteComputation
 * @section LICENCE
 *   Copyright (C) 1991 2023  EDF R&D                www.code-aster.org
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

#include "astercxx.h"

#include "aster_fort_calcul.h"
#include "aster_fort_superv.h"

#include "DataFields/FieldOnCellsBuilder.h"
#include "Discretization/Calcul.h"
#include "Discretization/DiscreteComputation.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Modeling/XfemModel.h"
#include "Utilities/Tools.h"

ElementaryVectorDisplacementRealPtr DiscreteComputation::getMechanicalNeumannForces(
    const ASTERDOUBLE time_curr, const ASTERDOUBLE time_step, const ASTERDOUBLE theta ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    auto elemVect = std::make_shared< ElementaryVectorDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_MECA" );
    elemVect->prepareCompute( calcul_option );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );

    FieldOnCellsRealPtr externVar = nullptr;
    if ( currMater && currMater->hasExternalStateVariable() ) {
        externVar = _phys_problem->getExternalStateVariables( time_curr );
    }

    ConstantFieldOnCellsRealPtr currBehav = nullptr;
    if ( currMater ) {
        currBehav = currMater->getBehaviourField();
    }

    auto isXfem = currModel->existsXfem();

    auto calcul = std::make_unique< Calcul >( calcul_option );

    auto impl = [&]( auto load, const ASTERINTEGER &load_i, const std::string &option,
                     const std::string &name, const std::string &param,
                     const FiniteElementDescriptorPtr FED,
                     std::vector< std::pair< std::string, DataFieldPtr > > field_in = {} ) {
        if ( load->hasLoadField( name ) ) {
            calcul->setOption( option );
            calcul->setFiniteElementDescriptor( FED );

            calcul->clearInputs();
            calcul->addTimeField( "PTEMPSR", time_curr, time_step, theta );
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            if ( currMater ) {
                calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
            }
            if ( currElemChara ) {
                calcul->addElementaryCharacteristicsField( currElemChara );
            }
            if ( isXfem ) {
                calcul->addXFEMField( currModel->getXfemModel() );
            }
            calcul->addInputField( param, load->getConstantLoadField( name ) );

            for ( auto &[param_in, field] : field_in ) {
                if ( field && field->exists() ) {
                    calcul->addInputField( param_in, field );
                }
            }

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ),
                                             load_i );
            }
        }
    };

    auto mecaLoadReal = listOfLoads->getMechanicalLoadsReal();
    for ( const auto &load : mecaLoadReal ) {
        auto load_FEDesc = load->getFiniteElementDescriptor();

        impl( load, iload, "CHAR_MECA_FORC_R", "FORNO", "PFORNOR", load_FEDesc );
        impl( load, iload, "CHAR_MECA_FR3D3D", "F3D3D", "PFR3D3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FRCO2D", "FCO2D", "PFRCO2D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FRCO3D", "FCO3D", "PFRCO3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FR2D3D", "F2D3D", "PFR2D3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FR1D3D", "F1D3D", "PFR1D3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FR2D2D", "F2D2D", "PFR2D2D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FR1D2D", "F1D2D", "PFR1D2D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FR1D1D", "F1D1D", "PFR1D1D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_PESA_R", "PESAN", "PPESANR", model_FEDesc );
        impl( load, iload, "CHAR_MECA_ROTA_R", "ROTAT", "PROTATR", model_FEDesc,
              {{"PCOMPOR", currBehav}} );
        impl( load, iload, "CHAR_MECA_EPSI_R", "EPSIN", "PEPSINR", model_FEDesc,
              {{"PCOMPOR", currBehav}} );
        impl( load, iload, "CHAR_MECA_FRELEC", "FELEC", "PFRELEC", model_FEDesc );
        impl( load, iload, "CHAR_MECA_PRES_R", "PRESS", "PPRESSR", model_FEDesc );
        impl( load, iload, "CHAR_MECA_ONDE", "ONDE", "PONDECR", model_FEDesc );

        iload++;
    }

    auto mecaLoadFunc = listOfLoads->getMechanicalLoadsFunction();
    for ( const auto &load : mecaLoadFunc ) {
        auto load_FEDesc = load->getFiniteElementDescriptor();

        impl( load, iload, "CHAR_MECA_FORC_F", "FORNO", "PFORNOF", load_FEDesc );
        impl( load, iload, "CHAR_MECA_FF3D3D", "F3D3D", "PFF3D3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FFCO2D", "FCO2D", "PFFCO2D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FFCO3D", "FCO3D", "PFFCO3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FF2D3D", "F2D3D", "PFF2D3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FF1D3D", "F1D3D", "PFF1D3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FF2D2D", "F2D2D", "PFF2D2D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FF1D2D", "F1D2D", "PFF1D2D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FF1D1D", "F1D1D", "PFF1D1D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_EPSI_F", "EPSIN", "PEPSINF", model_FEDesc,
              {{"PCOMPOR", currBehav}} );
        impl( load, iload, "CHAR_MECA_PRES_F", "PRESS", "PPRESSF", model_FEDesc );
        impl( load, iload, "ONDE_PLAN", "ONDPL", "PONDPLA", model_FEDesc,
              {{"PONDPLR", load->getConstantLoadField( "ONDPR" )}, {"PVARCPR", externVar}} );

        iload++;
    }

    elemVect->build();

    return elemVect;
};

/** @brief Compute internal forces, stress and internal state variables */
std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, FieldOnCellsRealPtr, FieldOnCellsRealPtr,
            FieldOnNodesRealPtr >
DiscreteComputation::getInternalForces( const FieldOnNodesRealPtr displ_prev,
                                        const FieldOnNodesRealPtr displ_step,
                                        const FieldOnCellsRealPtr stress,
                                        const FieldOnCellsRealPtr internVar,
                                        const ASTERDOUBLE &time_prev, const ASTERDOUBLE &time_step,
                                        const FieldOnCellsRealPtr &externVarPrev,
                                        const FieldOnCellsRealPtr &externVarCurr,
                                        const VectorString &groupOfCells ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currBehaviour = _phys_problem->getBehaviourProperty();

    // Select option to compute
    std::string option = "RAPH_MECA";

    // Prepare computing:
    CalculPtr calcul = createCalculForNonLinear( option, time_prev, time_prev + time_step,
                                                 externVarPrev, externVarCurr, groupOfCells );
    FiniteElementDescriptorPtr FEDesc = calcul->getFiniteElementDescriptor();

    // Set current physical state
    calcul->addInputField( "PDEPLMR", displ_prev );
    calcul->addInputField( "PDEPLPR", displ_step );
    calcul->addInputField( "PCONTMR", stress );
    calcul->addInputField( "PVARIMR", internVar );

    // Coded Material
    auto currCodedMater = _phys_problem->getCodedMaterial();
    calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

    // Provisoire: pour TANGENTE=VERIFICATION, nécessité de variables internes à chaque itération
    auto vari_iter = FieldOnCellsPtrBuilder< ASTERDOUBLE >( FEDesc, "ELGA", "VARI_R", currBehaviour,
                                                            currElemChara );
    calcul->addInputField( "PVARIMP", vari_iter );

    // Create output vector
    auto elemVect = std::make_shared< ElementaryVectorReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );
    elemVect->prepareCompute( option );

    // Create output fields
    auto stress_curr =
        FieldOnCellsPtrBuilder< ASTERDOUBLE >( FEDesc, "ELGA", "SIEF_R", currElemChara );
    FieldOnCellsLongPtr exitField = std::make_shared< FieldOnCellsLong >( FEDesc );
    auto vari_curr = FieldOnCellsPtrBuilder< ASTERDOUBLE >( FEDesc, "ELGA", "VARI_R", currBehaviour,
                                                            currElemChara );

    // Add output fields
    calcul->addOutputField( "PVARIPR", vari_curr );
    calcul->addOutputField( "PCONTPR", stress_curr );
    calcul->addOutputField( "PCODRET", exitField );

    // Add output elementary
    calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );

    // Compute and assemble vector
    FieldOnNodesRealPtr internalForces;
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) )
            elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ) );
        elemVect->build();
        internalForces = elemVect->assemble( _phys_problem->getDOFNumbering() );
    };

    std::string exitFieldName = ljust( exitField->getName(), 19 );
    ASTERINTEGER exitCode = 0;
    CALLO_GETERRORCODE( exitFieldName, &exitCode );
#ifdef ASTER_HAVE_MPI
    ASTERINTEGER exitCodeLocal = exitCode;
    AsterMPI::all_reduce( exitCodeLocal, exitCode, MPI_MAX );
#endif

    return std::make_tuple( exitField, exitCode, vari_curr, stress_curr, internalForces );
}

/** @brief Compute AFFE_CHAR_MECA DDL_IMPO */
ElementaryVectorDisplacementRealPtr
DiscreteComputation::getMechanicalImposedDualBC( const ASTERDOUBLE time_curr ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    auto elemVect = std::make_shared< ElementaryVectorDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_MECA" );
    elemVect->prepareCompute( calcul_option );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto listOfLoads = _phys_problem->getListOfLoads();

    auto calcul = std::make_unique< Calcul >( calcul_option );

    auto impl_disp = [&]( auto loads, bool real ) {
        std::string name;
        if ( real ) {
            calcul->setOption( "MECA_DDLI_R" );
            name = "PDDLIMR";
        } else {
            calcul->setOption( "MECA_DDLI_F" );
            name = "PDDLIMF";
        }
        for ( const auto &load : loads ) {
            auto load_FEDesc = load->getFiniteElementDescriptor();
            auto impo_field = load->getImposedField();
            if ( impo_field && impo_field->exists() && load_FEDesc ) {
                calcul->clearInputs();
                calcul->clearOutputs();
                calcul->setFiniteElementDescriptor( load_FEDesc );
                calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
                calcul->addInputField( name, impo_field );
                if ( !real ) {
                    calcul->addTimeField( "PTEMPSR", time_curr );
                }
                calcul->addOutputElementaryTerm( "PVECTUR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ),
                                                 iload );
                }
            }
            iload++;
        }
    };

    impl_disp( listOfLoads->getMechanicalLoadsReal(), true );
    impl_disp( listOfLoads->getMechanicalLoadsFunction(), false );

#ifdef ASTER_HAVE_MPI
    impl_disp( listOfLoads->getParallelMechanicalLoadsReal(), true );
    impl_disp( listOfLoads->getParallelMechanicalLoadsFunction(), false );
#endif

    auto impl_vite = [&]( auto loads, bool real ) {
        std::string name;
        if ( real ) {
            calcul->setOption( "CHAR_MECA_VFAC" );
            name = "PVITEFR";
        } else {
            calcul->setOption( "CHAR_MECA_VFAC_F" );
            name = "PVITEFF";
        }

        for ( const auto &load : loads ) {
            if ( load->hasLoadField( "VFACE" ) ) {
                calcul->clearInputs();
                calcul->clearOutputs();
                calcul->setFiniteElementDescriptor( currModel->getFiniteElementDescriptor() );
                calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
                auto currCodedMater = _phys_problem->getCodedMaterial();
                calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
                calcul->addInputField( name, load->getConstantLoadField( "VFACE" ) );
                if ( !real ) {
                    calcul->addTimeField( "PTEMPSR", time_curr );
                }
                calcul->addOutputElementaryTerm( "PVECTUR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ),
                                                 iload );
                }
            }
            iload++;
        }
    };

    impl_vite( listOfLoads->getMechanicalLoadsReal(), true );
    impl_vite( listOfLoads->getMechanicalLoadsFunction(), false );

    elemVect->build();

    return elemVect;
}

FieldOnNodesRealPtr DiscreteComputation::getContactForces(
    const MeshCoordinatesFieldPtr geom, const FieldOnNodesRealPtr displ_prev,
    const FieldOnNodesRealPtr displ_step, const ASTERDOUBLE &time_prev,
    const ASTERDOUBLE &time_step, const FieldOnCellsRealPtr data,
    const FieldOnNodesRealPtr coef_cont, const FieldOnNodesRealPtr coef_frot ) const {
    // Select option for matrix
    std::string option = "CHAR_MECA_CONT";

    auto [Fed_Slave, Fed_pair] = _phys_problem->getListOfLoads()->getContactLoadDescriptor();

    // Prepare computing
    CalculPtr calcul = std::make_unique< Calcul >( option );
    calcul->setFiniteElementDescriptor( Fed_pair );

    // Set input field
    calcul->addInputField( "PGEOMER", _phys_problem->getMesh()->getCoordinates() );
    calcul->addInputField( "PGEOMCR", geom );
    calcul->addInputField( "PDEPL_M", displ_prev );
    calcul->addInputField( "PDEPL_P", displ_step );
    calcul->addInputField( "PCONFR", data );
    calcul->addInputField( "PCCONTR", coef_cont );
    calcul->addInputField( "PCFROTR", coef_frot );

    // Add time fields
    calcul->addTimeField( "PINSTMR", time_prev );
    calcul->addTimeField( "PINSTPR", time_prev + time_step );

    // Create output vector
    auto elemVect = std::make_shared< ElementaryVectorDisplacementReal >();
    elemVect->prepareCompute( option );

    // Add output elementary
    calcul->addOutputElementaryTerm( "PVECTCR", std::make_shared< ElementaryTermReal >() );
    calcul->addOutputElementaryTerm( "PVECTFR", std::make_shared< ElementaryTermReal >() );

    // Computation
    calcul->compute();
    if ( calcul->hasOutputElementaryTerm( "PVECTCR" ) ) {
        elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTCR" ) );
    }
    if ( calcul->hasOutputElementaryTerm( "PVECTFR" ) ) {
        elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTFR" ) );
    }
    elemVect->build();

    return elemVect->assemble( _phys_problem->getDOFNumbering() );
}
