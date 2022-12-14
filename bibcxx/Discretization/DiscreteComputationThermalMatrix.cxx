/**
 * @file DiscreteComputation.cxx
 * @brief Implementation of class DiscreteComputation
 * @section LICENCE
 *   Copyright (C) 1991 2022  EDF R&D                www.code-aster.org
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

#include "aster_fort_calcul.h"
#include "aster_fort_superv.h"
#include "astercxx.h"

#include "Discretization/Calcul.h"
#include "Discretization/DiscreteComputation.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Loads/ThermalLoad.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Modeling/XfemModel.h"
#include "Utilities/Tools.h"

ElementaryMatrixTemperatureRealPtr DiscreteComputation::getLinearConductivityMatrix(
    const ASTERDOUBLE time, const ASTERINTEGER &modeFourier, const VectorString &groupOfCells,
    const bool &with_dual ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    const std::string option( "RIGI_THER" );

    auto elemMatr = std::make_shared< ElementaryMatrixTemperatureReal >( _phys_problem );
    elemMatr->prepareCompute( option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Prepare computing
    CalculPtr calcul = std::make_shared< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

        if ( currMater->hasExternalStateVariable() ) {
            calcul->addInputField( "PVARCPR", _phys_problem->getExternalStateVariables( time ) );
        }
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
    }

    calcul->addFourierModeField( modeFourier );
    calcul->addTimeField( "PTEMPSR", time, 1.0, 1.0 );

    if ( currModel->existsXfem() ) {
        XfemModelPtr currXfemModel = currModel->getXfemModel();
        calcul->addXFEMField( currXfemModel );
    }

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for rigidity
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
    };

    if ( with_dual ) {
        DiscreteComputation::baseDualLinearConductivityMatrix( calcul, elemMatr );
        DiscreteComputation::baseExchangeThermalMatrix( calcul, elemMatr, time );
    }

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixTemperatureRealPtr
DiscreteComputation::getLinearCapacityMatrix( const ASTERDOUBLE time,
                                              const VectorString &groupOfCells ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    const std::string option( "MASS_THER" );

    auto elemMatr = std::make_shared< ElementaryMatrixTemperatureReal >( _phys_problem );
    elemMatr->prepareCompute( option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Prepare computing
    auto calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    // Set to -1 because not used.
    calcul->addTimeField( "PTEMPSR", time, 1.0, -1.0 );

    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

        if ( currMater->hasExternalStateVariable() ) {
            calcul->addInputField( "PVARCPR", _phys_problem->getExternalStateVariables( time ) );
        }
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
    }

    if ( currModel->existsXfem() ) {
        XfemModelPtr currXfemModel = currModel->getXfemModel();
        calcul->addXFEMField( currXfemModel );
    }

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for mass
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
    };

    elemMatr->build();
    return elemMatr;
};

void DiscreteComputation::baseDualLinearConductivityMatrix(
    CalculPtr &calcul, ElementaryMatrixTemperatureRealPtr &elemMatr ) const {

    // Prepare loads
    const auto listOfLoads = _phys_problem->getListOfLoads();

    // Select option
    calcul->setOption( "THER_DDLM_R" );

    auto impl = [calcul, elemMatr]( auto loads ) {
        for ( const auto &load : loads ) {
            auto FEDesc = load->getFiniteElementDescriptor();
            auto field = load->getMultiplicativeField();
            if ( field && field->exists() && FEDesc ) {
                calcul->clearInputs();
                calcul->clearOutputs();
                calcul->setFiniteElementDescriptor( FEDesc );
                calcul->addInputField( "PDDLMUR", field );
                calcul->addOutputElementaryTerm( "PMATTTR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                    elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
                }
            }
        }
    };

    impl( listOfLoads->getThermalLoadsReal() );
    impl( listOfLoads->getThermalLoadsFunction() );

#ifdef ASTER_HAVE_MPI
    impl( listOfLoads->getParallelThermalLoadsReal() );
    impl( listOfLoads->getParallelThermalLoadsFunction() );
#endif
};

ElementaryMatrixTemperatureRealPtr DiscreteComputation::getDualLinearConductivityMatrix() const {
    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    const std::string option( "THER_DDLM_R" );

    auto elemMatr = std::make_shared< ElementaryMatrixTemperatureReal >( _phys_problem );
    elemMatr->prepareCompute( option );

    // Prepare computing
    CalculPtr calcul = std::make_unique< Calcul >( option );

    // Compute elementary matrices
    DiscreteComputation::baseDualLinearConductivityMatrix( calcul, elemMatr );

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixTemperatureRealPtr
DiscreteComputation::getExchangeThermalMatrix( const ASTERDOUBLE &time ) const {
    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    const std::string option( "RIGI_THER" );

    auto elemMatr = std::make_shared< ElementaryMatrixTemperatureReal >( _phys_problem );
    elemMatr->prepareCompute( option );

    // Prepare computing
    CalculPtr calcul = std::make_unique< Calcul >( option );

    // Compute elementary matrices
    DiscreteComputation::baseExchangeThermalMatrix( calcul, elemMatr, time );

    elemMatr->build();
    return elemMatr;
};

void DiscreteComputation::baseExchangeThermalMatrix( CalculPtr &calcul,
                                                     ElementaryMatrixTemperatureRealPtr &elemMatr,
                                                     const ASTERDOUBLE &time ) const {

    // Prepare loads
    const auto &_listOfLoads = _phys_problem->getListOfLoads();
    auto isXfem = _phys_problem->getModel()->existsXfem();

    auto therLoadReal = _listOfLoads->getThermalLoadsReal();
    for ( const auto &load : therLoadReal ) {
        auto FEDesc = _phys_problem->getModel()->getFiniteElementDescriptor();
        auto load_FEDesc = load->getFiniteElementDescriptor();

        if ( load->hasLoadResult() ) {
            std::string evol_char_name = load->getLoadResultName();
            FieldOnCellsRealPtr evol_exchange_field =
                std::make_shared< FieldOnCellsReal >( FEDesc );
            std::string para( "COEF_H" );
            std::string access_var( "INST" );
            std::string base( "G" );
            std::string extr_right( "EXCLU" );
            std::string extr_left( "EXCLU" );
            ASTERINTEGER iret = 100;
            ASTERINTEGER stop = 0;

            CALLO_RSINCH( evol_char_name, para, access_var, &time, evol_exchange_field->getName(),
                          extr_right, extr_left, &stop, base, &iret );

            if ( iret >= 2 ) {
                AS_ABORT( "Cannot find COEF_H in EVOL_CHAR " + evol_char_name + " at time " +
                          std::to_string( time ) );
            }

            calcul->setOption( "RIGI_THER_COEH_R" );
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( FEDesc );
            calcul->addInputField( "PGEOMER",
                                   _phys_problem->getModel()->getMesh()->getCoordinates() );
            calcul->addInputField( "PCOEFHR", evol_exchange_field );
            calcul->addTimeField( "PTEMPSR", time, -1.0, 1.0 );
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
            }
        }

        if ( load->hasLoadField( "COEFH" ) ) {
            auto exchange_field = load->getConstantLoadField( "COEFH" );
            calcul->setOption( "RIGI_THER_COEH_R" );
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( FEDesc );
            calcul->addInputField( "PGEOMER",
                                   _phys_problem->getModel()->getMesh()->getCoordinates() );
            calcul->addInputField( "PCOEFHR", exchange_field );
            calcul->addTimeField( "PTEMPSR", time, -1.0, 1.0 );
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
            }
        }

        if ( load->hasLoadField( "HECHP" ) ) {
            auto wall_exchange_field = load->getConstantLoadField( "HECHP" );
            calcul->setOption( "RIGI_THER_PARO_R" );
            calcul->clearInputs();
            calcul->clearOutputs();
            if ( isXfem ) {
                XfemModelPtr currXfemModel = _phys_problem->getModel()->getXfemModel();
                calcul->addXFEMField( currXfemModel );
                calcul->setFiniteElementDescriptor( FEDesc );
            } else {
                calcul->setFiniteElementDescriptor( load_FEDesc );
            }
            calcul->addInputField( "PGEOMER",
                                   _phys_problem->getModel()->getMesh()->getCoordinates() );
            calcul->addInputField( "PHECHPR", wall_exchange_field );
            calcul->addTimeField( "PTEMPSR", time, -1.0, 1.0 );

            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
            }
        }
    }

    auto therLoadFunc = _listOfLoads->getThermalLoadsFunction();
    for ( const auto &load : therLoadFunc ) {
        auto FEDesc = _phys_problem->getModel()->getFiniteElementDescriptor();
        auto load_FEDesc = load->getFiniteElementDescriptor();

        if ( load->hasLoadField( "COEFH" ) ) {
            auto exchange_field = load->getConstantLoadField( "COEFH" );
            calcul->setOption( "RIGI_THER_COEH_F" );
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( FEDesc );
            calcul->addInputField( "PGEOMER",
                                   _phys_problem->getModel()->getMesh()->getCoordinates() );
            calcul->addInputField( "PCOEFHF", exchange_field );
            calcul->addTimeField( "PTEMPSR", time, -1.0, 1.0 );
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
            }
        }

        if ( load->hasLoadField( "HECHP" ) ) {
            auto wall_exchange_field = load->getConstantLoadField( "HECHP" );
            calcul->setOption( "RIGI_THER_PARO_F" );
            calcul->clearInputs();
            calcul->clearOutputs();
            if ( isXfem ) {
                XfemModelPtr currXfemModel = _phys_problem->getModel()->getXfemModel();
                calcul->addXFEMField( currXfemModel );
                calcul->setFiniteElementDescriptor( FEDesc );
            } else {
                calcul->setFiniteElementDescriptor( load_FEDesc );
            }
            calcul->addInputField( "PGEOMER",
                                   _phys_problem->getModel()->getMesh()->getCoordinates() );
            calcul->addInputField( "PHECHPF", wall_exchange_field );
            calcul->addTimeField( "PTEMPSR", time, -1.0, 1.0 );
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
            }
        }
    }

    return;
}
