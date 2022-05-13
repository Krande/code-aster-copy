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
    const ASTERDOUBLE &time, const ASTERDOUBLE &delta_time, const ASTERINTEGER &modeFourier,
    const VectorString &groupOfCells, const FieldOnCellsRealPtr _externVarField ) {

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
    CalculPtr _calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        _calcul->setModel( currModel );
    } else {
        _calcul->setGroupsOfCells( currModel, groupOfCells );
    }
    elemMatr->prepareCompute( option );

    // Add input fields
    _calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    if ( currMater ) {
        _calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
    }
    if ( _externVarField ) {
        _calcul->addInputField( "PVARCPR", _externVarField );
    }
    if ( currElemChara ) {
        _calcul->addElementaryCharacteristicsField( currElemChara );
    }

    _calcul->addFourierModeField( modeFourier );
    _calcul->addTimeField( time, delta_time, 1., 0., 0., 0. );

    if ( currModel->existsXfem() ) {
        XfemModelPtr currXfemModel = currModel->getXfemModel();
        _calcul->addXFEMField( currXfemModel );
    }

    // Add output elementary terms
    _calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for rigidity
    if ( currModel->existsFiniteElement() ) {
        _calcul->compute();
        std::cout << _calcul->hasOutputElementaryTerm( "PMATTTR" ) << std::endl;
        if ( _calcul->hasOutputElementaryTerm( "PMATTTR" ) )
            elemMatr->addElementaryTerm( _calcul->getOutputElementaryTerm( "PMATTTR" ) );
    };

    // Compute elementary matrices for dual BC
    DiscreteComputation::baseDualThermalMatrix( _calcul, elemMatr );

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixTemperatureRealPtr
DiscreteComputation::linearCapacityMatrix( const ASTERDOUBLE &time,
                                           const VectorString &groupOfCells,
                                           const FieldOnCellsRealPtr _externVarField ) {
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

    // // Check super-element
    // if ( currModel->existsSuperElement() ) {
    //     std::string modelName = ljust( currModel->getName(), 8 );
    //     CALLO_CHECKSUPERELEMENT( option, modelName );
    // }

    // Prepare computing
    auto _calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        _calcul->setModel( currModel );
    } else {
        _calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    elemMatr->prepareCompute( option );

    // Add input fields
    _calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    if ( currMater ) {
        _calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
    }
    if ( _externVarField ) {
        _calcul->addInputField( "PVARCPR", _externVarField );
    }
    if ( currElemChara ) {
        _calcul->addElementaryCharacteristicsField( currElemChara );
    }

    if ( currModel->existsXfem() ) {
        XfemModelPtr currXfemModel = currModel->getXfemModel();
        _calcul->addXFEMField( currXfemModel );
    }

    // Add output elementary terms
    _calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for mass
    if ( currModel->existsFiniteElement() ) {
        _calcul->compute();
        if ( _calcul->hasOutputElementaryTerm( "PMATTTR" ) )
            elemMatr->addElementaryTerm( _calcul->getOutputElementaryTerm( "PMATTTR" ) );
    };

    elemMatr->build();
    return elemMatr;
};

void DiscreteComputation::baseDualThermalMatrix( CalculPtr &calcul,
                                                 ElementaryMatrixTemperatureRealPtr &elemMatr ) {

    // Prepare loads
    const auto &_listOfLoads = _phys_problem->getListOfLoads();

    // Select option
    calcul->setOption( "THER_DDLM_R" );
    calcul->clearInputs();
    calcul->clearOutputs();

    auto therLoadReal = _listOfLoads->getThermalLoadsReal();
    for ( const auto &curIter : therLoadReal ) {
        auto FEDesc = curIter->getFiniteElementDescriptor();
        auto field = curIter->getThermalLoadDescription()->getMultiplicativeField();
        if ( field && field->exists() && FEDesc ) {
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( FEDesc );
            calcul->addInputField( "PDDLMUR", field );
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATTTR" ) );
            }
        }
    }

    auto therLoadFunc = _listOfLoads->getThermalLoadsFunction();
    for ( const auto &curIter : therLoadFunc ) {
        auto FEDesc = curIter->getFiniteElementDescriptor();
        auto field = curIter->getThermalLoadDescription()->getMultiplicativeField();
        if ( field && field->exists() && FEDesc ) {
            auto resuElem = std::make_shared< ElementaryTermReal >();
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( FEDesc );
            calcul->addInputField( "PDDLMUR", field );
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATTTR" ) );
            }
        }
    }

    // #ifdef ASTER_HAVE_MPI
    //     auto mecaParaLoadReal = _listOfLoads->getParallelMechanicalLoadsReal();
    //     for ( const auto &curIter : mecaParaLoadReal ) {
    //         auto FEDesc = curIter->getFiniteElementDescriptor();
    //         auto field = curIter->getMultiplicativeField();
    //         if ( field && field->exists() && FEDesc ) {
    //             auto resuElem = std::make_shared< ElementaryTermReal >();
    //             calcul->clearInputs();
    //             calcul->clearOutputs();
    //             calcul->setFiniteElementDescriptor( FEDesc );
    //             calcul->addInputField( "PDDLMUR", field );
    //             calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal
    //             >() ); calcul->compute(); if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
    //                elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATTTR" ) );
    //            }
    //         }
    //     }

    //     auto mecaParaLoadFunc = _listOfLoads->getParallelMechanicalLoadsFunction();
    //     for ( const auto &curIter : mecaParaLoadFunc ) {
    //         auto FEDesc = curIter->getFiniteElementDescriptor();
    //         auto field = curIter->getMultiplicativeField();
    //         if ( field && field->exists() && FEDesc ) {
    //             auto resuElem = std::make_shared< ElementaryTermReal >();
    //             calcul->clearInputs();
    //             calcul->clearOutputs();
    //             calcul->setFiniteElementDescriptor( FEDesc );
    //             calcul->addInputField( "PDDLMUR", field );
    //             calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal
    //             >() ); calcul->compute();
    //            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
    //             elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATTTR" ) );
    //           }
    //         }
    //     }

    // #endif
};
