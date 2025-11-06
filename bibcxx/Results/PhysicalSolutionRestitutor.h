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

#pragma once

#include "astercxx.h"

#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"
#include "Results/GeneralizedResult.h"
#include "Results/ModeResult.h"
#include "Results/Result.h"
#include "Studies/PhysicalProblem.h"

#include <cblas.h>
#include <immintrin.h>
#include <omp.h>

class PhysicalSolutionRestitutor {
  private:
    VectorReal _gene_disp;
    VectorReal _gene_velo;
    VectorReal _gene_acce;
    ModeResultPtr _mode_resu;
    TransientGeneralizedResultPtr _gene_resu;
    ASTERINTEGER _ar;
    ASTERINTEGER _nbatch;

  public:
    // --- Constructeur ---
    PhysicalSolutionRestitutor( const TransientGeneralizedResultPtr &resgen,
                                const ASTERINTEGER &nbatch, const ASTERINTEGER &ar )
        : _gene_resu( resgen ), _nbatch( nbatch ), _ar( ar ) {
        if ( !_gene_resu ) {
            throw std::runtime_error( "PhysicalSolutionRestitutor: resgen is null." );
        }

        // Retrieve the ModeResultPtr directly from resgen
        ModeResultPtr modPtr =
            _gene_resu->getGeneralizedDOFNumbering()->getModalBasisFromModeResult();
        if ( !modPtr ) {
            throw std::runtime_error( "PhysicalSolutionRestitutor: modPtr (mode result) is null." );
        }

        _mode_resu = modPtr; // store it internally

        // Immediate loading of concatenated vectors
        _gene_disp = _gene_resu->getDisplacementValues();
        _gene_acce = _gene_resu->getAccelerationValues();
        _gene_velo = _gene_resu->getVelocityValues();

        if ( _gene_disp.empty() ) {
            throw std::runtime_error(
                "PhysicalSolutionRestitutor: displacement coefficients are empty." );
        }
    }

    // --- Accesseurs utiles ---
    const VectorReal &getDisplacementCoeffs() const noexcept { return _gene_disp; }
    const VectorReal &getAccelerationCoeffs() const noexcept { return _gene_acce; }
    const VectorReal &getVelocityCoeffs() const noexcept { return _gene_velo; }

    const ModeResultPtr &getModeResults() const noexcept { return _mode_resu; }
    const TransientGeneralizedResultPtr &getGeneralizedResults() const noexcept {
        return _gene_resu;
    }

    struct FieldIndexInfo {
        VectorLong start_indices;
        VectorLong ncomps;
    };

    inline void matVecBatch( const double *A_ptr, const double *X_ptr, double *Y_ptr, size_t nVals,
                             size_t nModes, size_t m ) {
        cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, nVals, m, nModes, 1.0, A_ptr, nVals,
                     X_ptr, nModes, 0.0, Y_ptr, nVals );
    }

    static inline void updateMaxAbsValues( const double *vals_in, double *vals_max,
                                           const size_t n ) {
        for ( size_t i = 0; i < n; ++i ) {
            const double vin = std::abs( vals_in[i] );
            vals_max[i] = std::max( vals_max[i], vin );
        }
    }

    static inline void updateMaxAbsValuesWithMoments( const double *vals_in, double *vals_max,
                                                      const FieldIndexInfo &field_info ) {
        const auto &indices = field_info.start_indices;
        const auto &ncmps = field_info.ncomps;
        const size_t nsp = indices.size();

        for ( size_t i = 0; i < nsp; ++i ) {
            const size_t idx = indices[i];
            const int ncmp = ncmps[i];

            vals_max[idx + 0] = std::max( vals_max[idx + 0], std::abs( vals_in[idx + 0] ) );
            vals_max[idx + 1] = std::max( vals_max[idx + 1], std::abs( vals_in[idx + 1] ) );
            vals_max[idx + 2] = std::max( vals_max[idx + 2], std::abs( vals_in[idx + 2] ) );

            // --- Moment components only if ncmp == 6 ---
            if ( ncmp == 6 ) {
                const double mx = vals_in[idx + 3];
                const double my = vals_in[idx + 4];
                const double mz = vals_in[idx + 5];
                const double norm_in = mx * mx + my * my + mz * mz;

                const double mx_old = vals_max[idx + 3];
                const double my_old = vals_max[idx + 4];
                const double mz_old = vals_max[idx + 5];
                const double norm_old = mx_old * mx_old + my_old * my_old + mz_old * mz_old;

                if ( norm_in > norm_old ) {
                    vals_max[idx + 3] = std::abs( mx );
                    vals_max[idx + 4] = std::abs( my );
                    vals_max[idx + 5] = std::abs( mz );
                }
            }
        }
    }

    template < typename FieldPtrType >
    static inline FieldIndexInfo computeFieldIndexInfo( const FieldPtrType &field ) {
        FieldIndexInfo info;

        if ( !field )
            throw std::runtime_error( "computeFieldIndexInfoGeneric: field is null." );

        const int ncomp_field = field->getNumberOfComponents();

        // Cas 1 : champ sur les nœuds
        if constexpr ( std::is_same_v< FieldPtrType, SimpleFieldOnNodesRealPtr > ) {
            const int nb_nodes = field->getNumberOfNodes();
            info.start_indices.reserve( nb_nodes );
            info.ncomps.reserve( nb_nodes );

            for ( ASTERINTEGER node = 0; node < nb_nodes; ++node ) {
                bool has_moment = true;

                // Détection des 3 ou 6 composantes
                for ( int j = 0; j < 6; ++j ) {
                    if ( !field->hasValue( node, j ) ) {
                        if ( j >= 3 )
                            has_moment = false;
                        break;
                    }
                }

                const int ncomp = has_moment ? 6 : 3;
                const size_t first_idx = static_cast< size_t >( node ) * ncomp_field;

                info.start_indices.push_back( first_idx );
                info.ncomps.push_back( ncomp );
            }
        }

        // Cas 2 : champ sur les cellules
        else if constexpr ( std::is_same_v< FieldPtrType, SimpleFieldOnCellsRealPtr > ) {
            const auto mesh = field->getMesh();
            if ( !mesh )
                throw std::runtime_error( "computeFieldIndexInfoGeneric: field mesh is null." );

            for ( const auto &cell : mesh->getCells() ) {
                ASTERINTEGER npt = field->getNumberOfPointsOfCell( cell );
                ASTERINTEGER nspt = field->getNumberOfSubPointsOfCell( cell );

                for ( ASTERINTEGER ipt = 0; ipt < npt; ++ipt ) {
                    for ( ASTERINTEGER ispt = 0; ispt < nspt; ++ispt ) {
                        bool has_moment = true;

                        for ( int j = 0; j < 6; ++j ) {
                            if ( !field->hasValue( cell, j, ipt, ispt ) ) {
                                if ( j >= 3 )
                                    has_moment = false;
                                break;
                            }
                        }

                        const int ncomp = has_moment ? 6 : 3;
                        const size_t first_idx = field->getPositionInArray( cell, 0, ipt, ispt );

                        info.start_indices.push_back( first_idx );
                        info.ncomps.push_back( ncomp );
                    }
                }
            }
        }

        else {
            throw std::runtime_error( "Unsupported field type in computeFieldIndexInfoGeneric." );
        }

        return info;
    }

    // computeTimeMaxOfModalFieldsOnNodes
    std::map< std::string, FieldOnNodesRealPtr > computeMaxForFieldsOnNodes();

    // computeTimeMaxOfModalFieldsOnCells
    std::map< std::string, FieldOnCellsRealPtr > computeMaxForFieldsOnCells();
};
