
/**
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

#include "DataFields/FieldConverter.h"
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

/** @brief Compute max value of EFGE_ELNO or EGRU_ELNO based on the maximum of the
 *         equivalent moment for piping studies*/
FieldOnCellsRealPtr
PostProcessing::computeMaxResultantForPipe( const ResultPtr &resu,
                                            const std::string &field_name ) const {

    std::string name = strip( field_name );

    if ( name != "EFGE_ELNO" && name != "EGRU_ELNO" ) {
        AS_ABORT( "Invalid field_name: " + name + ". Expected 'EFGE_ELNO' or 'EGRU_ELNO'." );
    }

    const auto &fieldNames = resu->getFieldsNames();
    if ( std::find( fieldNames.begin(), fieldNames.end(), name ) == fieldNames.end() ) {
        AS_ABORT( "Field " + name + " not found in resu->getFieldsNames()" );
    }

    VectorLong indexes = resu->getIndexesForFieldName( name );
    if ( indexes.empty() ) {
        AS_ABORT( "No index found for field " + name );
    }

    // Liste des champs (une par index)
    std::vector< SimpleFieldOnCellsRealPtr > input_fields;
    for ( auto idx : indexes ) {
        input_fields.push_back( toSimpleFieldOnCells( *resu->getFieldOnCellsReal( name, idx ) ) );
    }

    // Champ de sortie basé sur le premier champ
    auto field_out =
        std::make_shared< FieldOnCellsReal >( *resu->getFieldOnCellsReal( name, indexes[0] ) );
    field_out->updateValuePointers();
    auto field_out_s = toSimpleFieldOnCells( *field_out );

    // Boucle sur les mailles
    for ( const auto &cell : field_out_s->getMesh()->getCells() ) {
        ASTERINTEGER npt = field_out_s->getNumberOfPointsOfCell( cell );
        ASTERINTEGER nspt = field_out_s->getNumberOfSubPointsOfCell( cell );

        for ( ASTERINTEGER ipt = 0; ipt < npt; ++ipt ) {
            for ( ASTERINTEGER ispt = 0; ispt < nspt; ++ispt ) {

                double max_fx = -1.0;
                double max_fy = -1.0;
                double max_fz = -1.0;
                double max_moment = -1.0;

                double val_fx = std::numeric_limits< double >::quiet_NaN();
                double val_fy = std::numeric_limits< double >::quiet_NaN();
                double val_fz = std::numeric_limits< double >::quiet_NaN();

                std::array< double, 3 > moment_vec = { std::numeric_limits< double >::quiet_NaN(),
                                                       std::numeric_limits< double >::quiet_NaN(),
                                                       std::numeric_limits< double >::quiet_NaN() };

                for ( size_t i = 0; i < input_fields.size(); ++i ) {
                    double val[6];
                    bool has_moment = true;

                    for ( int j = 0; j < 6; ++j ) {
                        try {
                            val[j] = input_fields[i]->getValue( cell, j, ipt, ispt );
                        } catch ( ... ) {
                            val[j] = std::numeric_limits< double >::quiet_NaN();
                            if ( j >= 3 )
                                has_moment = false;
                        }
                    }

                    if ( !std::isnan( val[0] ) ) {
                        double fx = std::abs( val[0] );
                        if ( fx > max_fx ) {
                            max_fx = fx;
                            val_fx = fx;
                        }
                    }

                    if ( !std::isnan( val[1] ) ) {
                        double fy = std::abs( val[1] );
                        if ( fy > max_fy ) {
                            max_fy = fy;
                            val_fy = fy;
                        }
                    }

                    if ( !std::isnan( val[2] ) ) {
                        double fz = std::abs( val[2] );
                        if ( fz > max_fz ) {
                            max_fz = fz;
                            val_fz = fz;
                        }
                    }

                    if ( has_moment && !std::isnan( val[3] ) && !std::isnan( val[4] ) &&
                         !std::isnan( val[5] ) ) {
                        double mx = val[3], my = val[4], mz = val[5];
                        double norm = std::hypot( mx, my, mz );
                        if ( norm > max_moment ) {
                            max_moment = norm;
                            moment_vec = { std::abs( mx ), std::abs( my ), std::abs( mz ) };
                        }
                    }
                }

                if ( !std::isnan( val_fx ) || !std::isnan( val_fy ) || !std::isnan( val_fz ) ) {
                    field_out_s->setValue( cell, 0, ipt, ispt, val_fx ); // N
                    field_out_s->setValue( cell, 1, ipt, ispt, val_fy ); // VY
                    field_out_s->setValue( cell, 2, ipt, ispt, val_fz ); // VZ
                }
                if ( !std::isnan( moment_vec[0] ) ) {
                    field_out_s->setValue( cell, 3, ipt, ispt, moment_vec[0] ); // MT
                    field_out_s->setValue( cell, 4, ipt, ispt, moment_vec[1] ); // MFY
                    field_out_s->setValue( cell, 5, ipt, ispt, moment_vec[2] ); // MFZ
                }
            }
        }
    }

    return toFieldOnCells( *field_out_s, field_out->getDescription(), "EFGE_ELNO", "PEFFORR" );
}
