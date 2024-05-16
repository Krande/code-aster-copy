
/**
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

#include "PostProcessing/PostProcessing.h"

#include "DataFields/FieldOnCellsBuilder.h"
#include "Discretization/Calcul.h"
#include "Modeling/Model.h"

FieldOnCellsRealPtr PostProcessing::computeHydration( const FieldOnNodesRealPtr temp_prev,
                                                      const FieldOnNodesRealPtr temp_curr,
                                                      const ASTERDOUBLE time_prev,
                                                      const ASTERDOUBLE time_curr,
                                                      const FieldOnCellsRealPtr hydr_prev ) const {
    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    const std::string calcul_option( "HYDR_ELGA" );

    auto currModel = _phys_problem->getModel();
    auto currBehav = _phys_problem->getBehaviourProperty();
    auto currCodedMater = _phys_problem->getCodedMaterial();

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );

    // Add Input Field
    calcul->addInputField( "PCOMPOR", currBehav->getBehaviourField() );
    calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
    calcul->addTimeField( "PINSTR", time_curr, time_curr - time_prev, -1.0 );
    calcul->addInputField( "PTEMPMR", temp_prev );
    calcul->addInputField( "PTEMPPR", temp_curr );
    calcul->addInputField( "PHYDRMR", hydr_prev );

    // Create output fields
    auto hydr_curr = FieldOnCellsPtrBuilder< ASTERDOUBLE >( calcul->getFiniteElementDescriptor(),
                                                            "ELGA", "HYDR_R" );
    calcul->addOutputField( "PHYDRPR", hydr_curr );

    calcul->compute();

    return hydr_curr;
};

/** @brief Compute annealing */
FieldOnCellsRealPtr PostProcessing::computeAnnealing(
    const FieldOnCellsRealPtr internVar, const ASTERDOUBLE &time_prev, const ASTERDOUBLE &time_curr,
    const FieldOnCellsRealPtr &externVarPrev, const FieldOnCellsRealPtr &externVarCurr ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currBehaviour = _phys_problem->getBehaviourProperty();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Select option to compute
    std::string option = "REST_ECRO";

    // Prepare computing: the main object
    CalculPtr calcul = std::make_unique< Calcul >( option );
    calcul->setModel( currModel );

    // Add input fields
    calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
    calcul->addBehaviourField( currBehaviour );
    if ( currMater->hasExternalStateVariable() ) {
        if ( !externVarPrev ) {
            AS_ABORT( "External state variables vector for beginning of time step is missing" )
        }
        if ( !externVarCurr ) {
            AS_ABORT( "External state variables vector for end of time step is missing" )
        }
        calcul->addInputField( "PVARCMR", externVarPrev );
        calcul->addInputField( "PVARCPR", externVarCurr );
    }
    calcul->addTimeField( "PINSTMR", time_prev );
    calcul->addTimeField( "PINSTPR", time_curr );

    // Set current physical state
    calcul->addInputField( "PVARIMR", internVar );

    // Get Finite Element Descriptor
    FiniteElementDescriptorPtr FEDesc = calcul->getFiniteElementDescriptor();

    // Create output field
    auto vari_curr = FieldOnCellsPtrBuilder< ASTERDOUBLE >( FEDesc, "ELGA", "VARI_R", currBehaviour,
                                                            currElemChara );
    // Add output field
    calcul->addOutputField( "PVARIPR", vari_curr );

    // Compute
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
    };

    return vari_curr;
};
