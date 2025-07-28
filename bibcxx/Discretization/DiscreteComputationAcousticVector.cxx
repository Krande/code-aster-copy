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

#include "aster_fort_calcul.h"
#include "aster_fort_superv.h"

#include "Discretization/Calcul.h"
#include "Discretization/DiscreteComputation.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Modeling/XfemModel.h"
#include "Utilities/Tools.h"

std::variant< ElementaryVectorPressureComplexPtr, FieldOnNodesComplexPtr >
DiscreteComputation::getAcousticNeumannForces( const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isAcoustic() );

    auto elemVect =
        std::make_shared< ElementaryVectorPressureComplex >( _phys_problem->getModel() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_ACOU" );

    // Main parameters
    auto mesh = _phys_problem->getMesh();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currModel = _phys_problem->getModel();
    auto listOfLoads = _phys_problem->getListOfLoads();

    auto calcul = std::make_unique< Calcul >( calcul_option );

    auto acouLoadComplex = listOfLoads->getAcousticLoadsComplex();
    for ( const auto &load : acouLoadComplex ) {
        // Termes VITE
        if ( load->hasLoadField( "VFACE" ) ) {
            auto speed_field = load->getConstantLoadField( "VFACE" );

            calcul->setOption( "CHAR_ACOU_VFAC_C" );
            calcul->setFiniteElementDescriptor( currModel->getFiniteElementDescriptor() );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", mesh->getCoordinates() );
            calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
            calcul->addInputField( "PVITEFC", speed_field );

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTC",
                                             std::make_shared< ElementaryTermComplex >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTC" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermComplex( "PVECTTC" ),
                                             iload );
            }
        }

        iload++;
    }

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            raiseAsterError( "Not implemented" );
            return FieldOnNodesComplexPtr( nullptr );
        } else {
            auto vectAsse =
                std::make_shared< FieldOnNodesComplex >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

std::variant< ElementaryVectorPressureComplexPtr, FieldOnNodesComplexPtr >
DiscreteComputation::getAcousticVolumetricForces( const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isAcoustic() );

    auto elemVect =
        std::make_shared< ElementaryVectorPressureComplex >( _phys_problem->getModel() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_ACOU" );

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            raiseAsterError( "Not implemented" );
            return FieldOnNodesComplexPtr( nullptr );
        } else {
            auto vectAsse =
                std::make_shared< FieldOnNodesComplex >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

/** @brief Compute AFFE_CHAR_ACOU ACOU_IMPO */
std::variant< ElementaryVectorPressureComplexPtr, FieldOnNodesComplexPtr >
DiscreteComputation::getAcousticImposedDualBC( const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isAcoustic() );

    auto elemVect =
        std::make_shared< ElementaryVectorPressureComplex >( _phys_problem->getModel() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_ACOU" );

    // Main parameters
    auto mesh = _phys_problem->getMesh();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto listOfLoads = _phys_problem->getListOfLoads();

    auto calcul = std::make_unique< Calcul >( calcul_option );

    auto impl = [&]( auto loads, bool real ) {
        std::string name;
        if ( real ) {
            calcul->setOption( "ACOU_DDLI_C" );
            name = "PDDLIMC";
        } else {
            calcul->setOption( "ACOU_DDLI_F" );
            name = "PDDLIMF";
        }
        for ( const auto &load : loads ) {
            auto load_FEDesc = load->getFiniteElementDescriptor();
            auto impo_field = load->getImposedField();
            if ( impo_field && impo_field->exists() && load_FEDesc ) {
                calcul->clearInputs();
                calcul->clearOutputs();
                calcul->setFiniteElementDescriptor( load_FEDesc );
                calcul->addInputField( "PGEOMER", mesh->getCoordinates() );
                calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
                calcul->addInputField( name, impo_field );
                calcul->addOutputElementaryTerm( "PVECTTC",
                                                 std::make_shared< ElementaryTermComplex >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTTC" ) ) {
                    elemVect->addElementaryTerm(
                        calcul->getOutputElementaryTermComplex( "PVECTTC" ), iload );
                }
            }
            iload++;
        }
    };

    impl( listOfLoads->getAcousticLoadsComplex(), true );

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            raiseAsterError( "Not implemented" );
            return FieldOnNodesComplexPtr( nullptr );
        } else {
            auto vectAsse =
                std::make_shared< FieldOnNodesComplex >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};
