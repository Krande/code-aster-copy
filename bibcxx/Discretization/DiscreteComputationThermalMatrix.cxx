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

ElementaryMatrixTemperatureRealPtr DiscreteComputation::linearConductivityMatrix(
    const ASTERDOUBLE time_value, const ASTERDOUBLE time_delta, const ASTERDOUBLE time_theta,
    const ASTERINTEGER &modeFourier, const VectorString &groupOfCells,
    const FieldOnCellsRealPtr _externVarField ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemMatr = std::make_shared< ElementaryMatrixTemperatureReal >();

    // Get main parameters
    const std::string option( "RIGI_THER" );
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Check external state variables
    if ( currMater && currMater->hasExternalStateVariable() ) {
        if ( !_externVarField ) {
            AS_ABORT( "External state variables vector is missing" )
        }
    }

    // Set parameters of elementary matrix
    elemMatr->setModel( currModel );
    elemMatr->setMaterialField( currMater );
    elemMatr->setElementaryCharacteristics( currElemChara );

    // Prepare computing
    CalculPtr calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }
    elemMatr->prepareCompute( option );

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
    }
    if ( _externVarField ) {
        calcul->addInputField( "PVARCPR", _externVarField );
    }
    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
    }

    calcul->addFourierModeField( modeFourier );
    calcul->addTimeField( "PTEMPSR", time_value, time_delta, time_theta );

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
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATTTR" ) );
    };

    // Compute elementary matrices for dual BC
    DiscreteComputation::baseDualThermalMatrix( calcul, elemMatr );
    DiscreteComputation::baseExchangeThermalMatrix( calcul, elemMatr, time_value, time_delta,
                                                    time_theta );

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixTemperatureRealPtr DiscreteComputation::linearCapacityMatrix(
    const ASTERDOUBLE time_value, const ASTERDOUBLE time_delta, const ASTERDOUBLE time_theta,
    const VectorString &groupOfCells, const FieldOnCellsRealPtr _externVarField ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    auto elemMatr = std::make_shared< ElementaryMatrixTemperatureReal >();
    // Get main parameters
    const std::string option( "MASS_THER" );
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Check external state variables
    if ( currMater && currMater->hasExternalStateVariable() ) {
        if ( !_externVarField ) {
            AS_ABORT( "External state variables vector is missing" )
        }
    }

    // Set parameters of elementary matrix
    elemMatr->setModel( currModel );
    elemMatr->setMaterialField( currMater );
    elemMatr->setElementaryCharacteristics( currElemChara );

    // Prepare computing
    auto calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    elemMatr->prepareCompute( option );

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    calcul->addTimeField( "PTEMPSR", time_value, time_delta, time_theta );

    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
    }
    if ( _externVarField ) {
        calcul->addInputField( "PVARCPR", _externVarField );
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
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATTTR" ) );
    };

    elemMatr->build();
    return elemMatr;
};

void DiscreteComputation::baseDualThermalMatrix(
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
                    elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATTTR" ) );
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

void DiscreteComputation::baseExchangeThermalMatrix( CalculPtr &calcul,
                                                     ElementaryMatrixTemperatureRealPtr &elemMatr,
                                                     const ASTERDOUBLE time_value,
                                                     const ASTERDOUBLE time_delta,
                                                     const ASTERDOUBLE time_theta ) const {

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

            CALLO_RSINCH( evol_char_name, para, access_var, &time_value,
                          evol_exchange_field->getName(), extr_right, extr_left, &stop, base,
                          &iret );

            if ( iret >= 2 ) {
                AS_ABORT( "Cannot find COEF_H in EVOL_CHAR " + evol_char_name + " at time " +
                          std::to_string( time_value ) );
            }

            calcul->setOption( "RIGI_THER_COEH_R" );
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( FEDesc );
            calcul->addInputField( "PGEOMER",
                                   _phys_problem->getModel()->getMesh()->getCoordinates() );
            calcul->addInputField( "PCOEFHR", evol_exchange_field );
            calcul->addTimeField( "PTEMPSR", time_value, time_delta, time_theta );
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATTTR" ) );
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
            calcul->addTimeField( "PTEMPSR", time_value, time_delta, time_theta );
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATTTR" ) );
            }
        }

        if ( load->hasLoadField( "HECHP" ) ) {
            auto wall_exchange_field = load->getConstantLoadField( "HECHP" );
            calcul->setOption( "RIGI_THER_PARO_R" );
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->addInputField( "PGEOMER",
                                   _phys_problem->getModel()->getMesh()->getCoordinates() );
            calcul->addInputField( "PHECHPR", wall_exchange_field );
            calcul->addTimeField( "PTEMPSR", time_value, time_delta, time_theta );
            if ( isXfem ) {
                XfemModelPtr currXfemModel = _phys_problem->getModel()->getXfemModel();
                calcul->addXFEMField( currXfemModel );
                calcul->setFiniteElementDescriptor( FEDesc );
            } else {
                calcul->setFiniteElementDescriptor( load_FEDesc );
            }
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATTTR" ) );
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
            calcul->addTimeField( "PTEMPSR", time_value, time_delta, time_theta );
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATTTR" ) );
            }
        }

        if ( load->hasLoadField( "HECHP" ) ) {
            auto wall_exchange_field = load->getConstantLoadField( "HECHP" );
            calcul->setOption( "RIGI_THER_PARO_F" );
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->addInputField( "PGEOMER",
                                   _phys_problem->getModel()->getMesh()->getCoordinates() );
            calcul->addInputField( "PHECHPF", wall_exchange_field );
            calcul->addTimeField( "PTEMPSR", time_value, time_delta, time_theta );
            if ( isXfem ) {
                XfemModelPtr currXfemModel = _phys_problem->getModel()->getXfemModel();
                calcul->addXFEMField( currXfemModel );
                calcul->setFiniteElementDescriptor( FEDesc );
            } else {
                calcul->setFiniteElementDescriptor( load_FEDesc );
            }
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATTTR" ) );
            }
        }
    }

    return;
}
