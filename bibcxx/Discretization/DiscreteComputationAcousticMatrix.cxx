/**
 * @file DiscreteComputation.cxx
 * @brief Implementation of class DiscreteComputation
 * @section LICENCE
 *   Copyright (C) 1991 2025  EDF R&D                www.code-aster.org
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

#include "Discretization/Calcul.h"
#include "Discretization/DiscreteComputation.h"

void DiscreteComputation::baseDualAcousticMatrix(
    CalculPtr &calcul, ElementaryMatrixPressureComplexPtr &elemMatr ) const {

    // Prepare loads
    const auto listOfLoads = _phys_problem->getListOfLoads();

    // Select option
    calcul->setOption( "ACOU_DDLM_C" );

    auto impl = [calcul, elemMatr]( auto loads ) {
        for ( const auto &load : loads ) {
            auto FEDesc = load->getFiniteElementDescriptor();
            auto field = load->getMultiplicativeField();
            if ( field && field->exists() && FEDesc && FEDesc->exists() ) {
                calcul->clearInputs();
                calcul->clearOutputs();
                calcul->setFiniteElementDescriptor( FEDesc );
                calcul->addInputField( "PDDLMUC", field );
                calcul->addOutputElementaryTerm( "PMATTTC",
                                                 std::make_shared< ElementaryTermComplex >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PMATTTC" ) ) {
                    elemMatr->addElementaryTerm(
                        calcul->getOutputElementaryTermComplex( "PMATTTC" ) );
                }
            }
        }
    };

    impl( listOfLoads->getAcousticLoadsComplex() );
};

ElementaryMatrixPressureComplexPtr
DiscreteComputation::getLinearMobilityMatrix( const VectorString &groupOfCells,
                                              const bool &with_dual ) const {
    AS_ASSERT( _phys_problem->getModel()->isAcoustic() );

    const std::string option( "RIGI_ACOU" );

    auto elemMatr =
        std::make_shared< ElementaryMatrixPressureComplex >( _phys_problem->getModel(), option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currCodedMater = _phys_problem->getCodedMaterial();

    // Prepare computing
    auto calcul = std::make_shared< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PMATTTC", std::make_shared< ElementaryTermComplex >() );

    // Compute elementary matrices for mass
    calcul->compute();
    if ( calcul->hasOutputElementaryTerm( "PMATTTC" ) )
        elemMatr->addElementaryTerm( calcul->getOutputElementaryTermComplex( "PMATTTC" ) );

    if ( with_dual ) {
        DiscreteComputation::baseDualAcousticMatrix( calcul, elemMatr );
    }

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixPressureComplexPtr DiscreteComputation::getDualLinearMobilityMatrix() const {
    AS_ASSERT( _phys_problem->getModel()->isAcoustic() );

    const std::string option( "ACOU_DDLM_C" );

    auto elemMatr =
        std::make_shared< ElementaryMatrixPressureComplex >( _phys_problem->getModel(), option );

    // Prepare computing
    auto calcul = std::make_shared< Calcul >( option );
    calcul->setModel( _phys_problem->getModel() );

    DiscreteComputation::baseDualAcousticMatrix( calcul, elemMatr );

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixPressureComplexPtr
DiscreteComputation::getCompressibilityMatrix( const VectorString &groupOfCells ) const {

    AS_ASSERT( _phys_problem->getModel()->isAcoustic() );

    const std::string option( "MASS_ACOU" );

    auto elemMatr =
        std::make_shared< ElementaryMatrixPressureComplex >( _phys_problem->getModel(), option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currCodedMater = _phys_problem->getCodedMaterial();

    // Prepare computing
    auto calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PMATTTC", std::make_shared< ElementaryTermComplex >() );

    // Compute elementary matrices for mass
    calcul->compute();
    if ( calcul->hasOutputElementaryTerm( "PMATTTC" ) )
        elemMatr->addElementaryTerm( calcul->getOutputElementaryTermComplex( "PMATTTC" ) );

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixPressureComplexPtr
DiscreteComputation::getImpedanceMatrix( const ASTERINTEGER &onde_flui ) const {

    AS_ASSERT( _phys_problem->getModel()->isAcoustic() );

    const std::string option( "AMOR_ACOU" );

    auto elemMatr =
        std::make_shared< ElementaryMatrixPressureComplex >( _phys_problem->getModel(), option );

    // Get main parameters
    auto currMesh = _phys_problem->getMesh();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currModel = _phys_problem->getModel();

    auto calcul = std::make_unique< Calcul >( option );
    calcul->setModel( currModel );

    // Add input fields
    calcul->addInputField( "PGEOMER", currMesh->getCoordinates() );
    calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

    auto carteFluid = createWaveTypeFluidField( onde_flui );
    calcul->addInputField( "PWATFLAC", carteFluid );

    calcul->addOutputElementaryTerm( "PMATTTC", std::make_shared< ElementaryTermComplex >() );
    calcul->compute();
    if ( calcul->hasOutputElementaryTerm( "PMATTTC" ) ) {
        elemMatr->addElementaryTerm( calcul->getOutputElementaryTermComplex( "PMATTTC" ) );
    }

    elemMatr->build();
    return elemMatr;
};
