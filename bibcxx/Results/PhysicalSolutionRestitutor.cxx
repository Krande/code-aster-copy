/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#include "astercxx.h"

#include "Results/PhysicalSolutionRestitutor.h"

std::map< std::string, FieldOnNodesRealPtr > PhysicalSolutionRestitutor::computeMaxForFieldsOnNodes() {
    const auto &resmod = getModeResults();
    const auto &resgen = getGeneralizedResults();
    if ( !resmod || !resgen )
        throw std::runtime_error( "computeMaxForFieldsOnNodes: resmod or resgen is null." );

    const size_t m = _nbatch; // batch size
    const size_t nSteps = resgen->getIndexes().size();

    const VectorReal &coeffs_disp = getDisplacementCoeffs();
    const VectorReal &coeffs_vite = getVelocityCoeffs();
    const VectorReal &coeffs_acce = getAccelerationCoeffs();

    std::map< std::string, FieldOnNodesRealPtr > results;

    const std::vector< std::string > availableFields = resmod->getFieldsNames();
    auto hasField = [&]( const std::string &name ) {
        return std::find( availableFields.begin(), availableFields.end(), name ) !=
               availableFields.end();
    };

    // -------------------------------------------------------
    // Fonction générique pour les champs normaux (DEPL, VITE, ACCE)
    // -------------------------------------------------------
    auto computeGenericField = [&]( const std::string &name, const VectorReal &coeffs ) {
        VectorLong mode_idxs = resmod->getIndexesForFieldName( name );
        const size_t nModes = mode_idxs.size();
        if ( nModes == 0 )
            return FieldOnNodesRealPtr();

        std::vector< const double * > A_ptrs;
        A_ptrs.reserve( nModes );
        for ( auto idx : mode_idxs )
            A_ptrs.push_back( resmod->getFieldOnNodesReal( name, idx )->getValues()->getDataPtr() );

        auto field_ref = resmod->getFieldOnNodesReal( name, mode_idxs[0] );
        auto field_max = std::make_shared< FieldOnNodesReal >( *field_ref );

        double *max_ptr = field_max->getValues()->getDataPtr();
        const size_t nVals = field_max->getValues()->size();
        std::fill_n( max_ptr, nVals, 0.0 );

        std::vector< double > X_batch( m * nModes );
        std::vector< double > Y_batch( m * nVals );
        std::vector< double > A( nVals * nModes );
        for ( size_t j = 0; j < nModes; ++j )
            std::copy( A_ptrs[j], A_ptrs[j] + nVals, A.data() + j * nVals );

        const double *coeffs_ptr = coeffs.data();

        for ( size_t step = 0; step < nSteps; step += m ) {
            const size_t nBatch = std::min( m, nSteps - step );
            for ( size_t k = 0; k < nBatch; ++k )
                std::copy( coeffs_ptr + ( step + k ) * nModes,
                           coeffs_ptr + ( step + k + 1 ) * nModes, X_batch.data() + k * nModes );

            matVecBatch( A.data(), X_batch.data(), Y_batch.data(), nVals, nModes, nBatch );
            for ( size_t k = 0; k < nBatch; ++k )
                updateMaxAbsValues( Y_batch.data() + k * nVals, max_ptr, nVals );
        }

        return field_max;
    };

    // -------------------------------------------------------
    // Cas particulier : REAC_NODA avec SimpleFieldOnNodes
    // -------------------------------------------------------
    auto computeReacNodaField = [&]() {
        VectorLong mode_idxs = resmod->getIndexesForFieldName( "REAC_NODA" );
        const size_t nModes = mode_idxs.size();
        if ( nModes == 0 )
            return FieldOnNodesRealPtr();

        std::vector< const double * > A_ptrs;
        std::vector< SimpleFieldOnNodesRealPtr > fields;
        fields.reserve( nModes );

        for ( auto idx : mode_idxs ) {
            auto f = toSimpleFieldOnNodes( *resmod->getFieldOnNodesReal( "REAC_NODA", idx ) );
            A_ptrs.push_back( f->getValues()->getDataPtr() );
            fields.push_back( f );
        }

        auto field_ref =
            toSimpleFieldOnNodes( *resmod->getFieldOnNodesReal( "REAC_NODA", mode_idxs[0] ) );
        auto field_max =
            toSimpleFieldOnNodes( *resmod->getFieldOnNodesReal( "REAC_NODA", mode_idxs[0] ) );

        double *max_ptr = field_max->getValues()->getDataPtr();
        const size_t nVals = field_max->getValues()->size();
        std::fill_n( max_ptr, nVals, 0.0 );

        std::vector< double > X_batch( m * nModes );
        std::vector< double > Y_batch( m * nVals );
        std::vector< double > A( nVals * nModes );
        for ( size_t j = 0; j < nModes; ++j )
            std::copy( A_ptrs[j], A_ptrs[j] + nVals, A.data() + j * nVals );

        const double *coeffs_ptr = coeffs_disp.data();
        auto field_info = computeFieldIndexInfo( field_ref );

        for ( size_t step = 0; step < nSteps; step += m ) {
            const size_t nBatch = std::min( m, nSteps - step );
            for ( size_t k = 0; k < nBatch; ++k )
                std::copy( coeffs_ptr + ( step + k ) * nModes,
                           coeffs_ptr + ( step + k + 1 ) * nModes, X_batch.data() + k * nModes );

            matVecBatch( A.data(), X_batch.data(), Y_batch.data(), nVals, nModes, nBatch );
            for ( size_t k = 0; k < nBatch; ++k )
                updateMaxAbsValuesWithMoments( Y_batch.data() + k * nVals, max_ptr, field_info );
        }

        return toFieldOnNodes( *field_max );
    };

    if ( hasField( "DEPL" ) ) {
        results["DEPL"] = computeGenericField( "DEPL", coeffs_disp );
        results["VITE"] = computeGenericField( "DEPL", coeffs_vite );
        results["ACCE"] = computeGenericField( "DEPL", coeffs_acce );
    }
    if ( hasField( "REAC_NODA" ) )
        results["REAC_NODA"] =
            _ar ? computeReacNodaField() : computeGenericField( "REAC_NODA", coeffs_disp );

    // Optional: log missing fields
    for ( const std::string &f : { "DEPL", "REAC_NODA" } ) {
        if ( !hasField( f ) )
            std::cout << "[PhysicalSolutionRestitutor] Field not found in resmod: " << f << std::endl;
    }

    return results;
}

std::map< std::string, FieldOnCellsRealPtr > PhysicalSolutionRestitutor::computeMaxForFieldsOnCells() {
    const VectorReal &all_coeffs_disp = getDisplacementCoeffs();
    const ModeResultPtr &resmod = getModeResults();
    const TransientGeneralizedResultPtr &resgen = getGeneralizedResults();

    if ( !resmod || !resgen )
        throw std::runtime_error( "computeMaxForFieldsOnCells: resmod or resgen is null." );

    // --- Fields to process ---
    std::vector< std::string > field_names = { "EFGE_ELNO", "EGRU_ELNO" };

    std::map< std::string, FieldOnCellsRealPtr > results;

    for ( const auto &field_name : field_names ) {
        VectorLong mode_indexes = resmod->getIndexesForFieldName( field_name );
        const size_t nModes = mode_indexes.size();
        if ( nModes == 0 )
            continue; // skip unavailable field

        // --- Prepare pointers to modal field data ---
        std::vector< SimpleFieldOnCellsRealPtr >
            field_objects; // stocke les champs pour garder leur mémoire en vie
        std::vector< const double * > A_ptrs;

        for ( auto idx : mode_indexes ) {
            auto field_ptr =
                toSimpleFieldOnCells( *resmod->getFieldOnCellsReal( field_name, idx ) );
            field_objects.push_back( field_ptr ); // garde l'objet en vie
            A_ptrs.push_back( field_ptr->getValues()->getDataPtr() );
        }
        // --- Reference and output fields ---
        FieldOnCellsRealPtr field_ref = resmod->getFieldOnCellsReal( field_name, mode_indexes[0] );

        auto field_out = toSimpleFieldOnCells( *field_ref );
        auto field_max = toSimpleFieldOnCells( *field_ref );

        // --- Access data pointers once ---
        double *out_ptr = field_out->getValues()->getDataPtr();
        double *max_ptr = field_max->getValues()->getDataPtr();
        const size_t nVals = field_out->getValues()->size();

        std::memset( out_ptr, 0, nVals * sizeof( ASTERDOUBLE ) );
        std::memset( max_ptr, 0, nVals * sizeof( ASTERDOUBLE ) );

        // --- Time loop ---
        VectorReal coeff_disp( nModes );
        const VectorLong &gen_indexes = resgen->getIndexes();
        const size_t nSteps = gen_indexes.size();
        const double *all_coeffs_ptr = all_coeffs_disp.data();

        //

        const size_t m = _nbatch; // taille du batch
        std::vector< double > X_batch( m * nModes );
        std::vector< double > Y_batch( m * nVals );

        // Préparer A une seule fois (hors boucle)
        std::vector< double > A( nVals * nModes );
        for ( size_t j = 0; j < nModes; ++j ) {
            std::copy( A_ptrs[j], A_ptrs[j] + nVals, A.data() + j * nVals );
        }

        if ( _ar ) {
            auto field_info = computeFieldIndexInfo( field_out );

            for ( size_t step = 0; step < nSteps; step += m ) {
                const size_t nBatch = std::min( m, nSteps - step );

                // Charger les coefficients de batch
                for ( size_t k = 0; k < nBatch; ++k ) {
                    const double *coeff_step = all_coeffs_ptr + ( step + k ) * nModes;
                    std::copy( coeff_step, coeff_step + nModes, X_batch.data() + k * nModes );
                }

                // Produit vectorisé
                matVecBatch( A.data(), X_batch.data(), Y_batch.data(), nVals, nModes, nBatch );

                // Réduction max sur les nBatch résultats
                for ( size_t k = 0; k < nBatch; ++k ) {
                    const double *out_local = Y_batch.data() + k * nVals;
                    updateMaxAbsValuesWithMoments( out_local, max_ptr, field_info );
                }
            }
        } else {
            for ( size_t step = 0; step < nSteps; step += m ) {
                const size_t nBatch = std::min( m, nSteps - step );

                // Charger les coefficients
                for ( size_t k = 0; k < nBatch; ++k ) {
                    const double *coeff_step = all_coeffs_ptr + ( step + k ) * nModes;
                    std::copy( coeff_step, coeff_step + nModes, X_batch.data() + k * nModes );
                }

                // Produit vectorisé
                matVecBatch( A.data(), X_batch.data(), Y_batch.data(), nVals, nModes, nBatch );

                // Calcul du max
                for ( size_t k = 0; k < nBatch; ++k ) {
                    const double *out_ptr_local = Y_batch.data() + k * nVals;
                    updateMaxAbsValues( out_ptr_local, max_ptr, nVals );
                }
            }
        }

        results[field_name] =
            toFieldOnCells( *field_max, field_ref->getDescription(), "EFGE_ELNO", "PEFFORR" );
    }

    return results;
}