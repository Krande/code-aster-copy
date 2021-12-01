#ifndef ASTERMPI_H_
#define ASTERMPI_H_

/**
 * @file AsterMPI.h
 * @brief Fichier entete contenant des utilitaires de manipulation de containers STL en parallèle
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

/* person_in_charge: nicolas.sellenet at edf.fr */
#include <algorithm>
#include <numeric>
#include <set>

#include "aster_mpi.h"

#ifdef ASTER_HAVE_MPI

#include "astercxx.h"

#include "MemoryManager/JeveuxString.h"
#include "MemoryManager/JeveuxVector.h"

#endif /* ASTER_HAVE_MPI */

/** @brief Get MPI number of procs */
int getMPISize( aster_comm_t *comm = aster_get_current_comm() );

/** @brief Get MPI rank */
int getMPIRank( aster_comm_t *comm = aster_get_current_comm() );

#ifdef ASTER_HAVE_MPI

/* this code is inspired by the MPI class of the dolfin project : fenicsproject.org */

class AsterMPI {
  private:
    template < typename T > struct dependent_false : std::false_type {};
    template < typename T > static MPI_Datatype mpi_type() {
        static_assert( dependent_false< T >::value, "Unknown MPI type" );
        throw std::runtime_error( "Unknown MPI type" );
        return MPI_CHAR;
    }

  public:
    /// Gather values from each process (variable count per process)
    template < typename T >
    static void all_gather( const std::vector< T > &in_values,
                            std::vector< std::vector< T > > &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Gather values from each process (variable count per process)
    template < typename T >
    static void all_gather( const std::vector< T > &in_values, std::vector< T > &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Gather values from each process (variable count per process)
    template < typename T >
    static void all_gather( const std::vector< T > &in_values, JeveuxVector< T > &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Gather values, one primitive from each process (MPI_Allgather)
    template < typename T >
    static void all_gather( const T in_value, std::vector< T > &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Gather values, one primitive from each process (MPI_Allgather).
    /// Specialization for std::string
    static void all_gather( const std::string &in_values, VectorString &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Gather values, one primitive from each process (MPI_Allgather).
    /// Specialization for VectorString
    static void all_gather( const VectorString &in_values, VectorString &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Gather values, one primitive from each process (MPI_Allgather).
    /// Specialization for JeveuxString
    template < int length >
    static void all_gather( const std::vector< JeveuxString< length > > &in_values,
                            std::vector< JeveuxString< length > > &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// AllReduce values, one value from each process (MPI_AllReduce).
    template < typename T >
    static void all_reduce( const T in_value, T &out_value, MPI_Op op,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Broadcast one value from root
    template < typename T >
    static void bcast( T &value, int root, aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Broadcast a vector from root
    template < typename T >
    static void bcast( std::vector< T > &value, int root,
                       aster_comm_t *_commCurrent = aster_get_current_comm() );
};

//---------------------------------------------------------------------------
template <> inline MPI_Datatype AsterMPI::mpi_type< char >() { return MPI_CHAR; }
template <> inline MPI_Datatype AsterMPI::mpi_type< float >() { return MPI_FLOAT; }
template <> inline MPI_Datatype AsterMPI::mpi_type< double >() { return MPI_DOUBLE; }
template <> inline MPI_Datatype AsterMPI::mpi_type< short int >() { return MPI_SHORT; }
template <> inline MPI_Datatype AsterMPI::mpi_type< int >() { return MPI_INT; }
template <> inline MPI_Datatype AsterMPI::mpi_type< long int >() { return MPI_LONG; }
template <> inline MPI_Datatype AsterMPI::mpi_type< unsigned int >() { return MPI_UNSIGNED; }
template <> inline MPI_Datatype AsterMPI::mpi_type< unsigned long int >() {
    return MPI_UNSIGNED_LONG;
}
template <> inline MPI_Datatype AsterMPI::mpi_type< long long >() { return MPI_LONG_LONG; }
template <> inline MPI_Datatype AsterMPI::mpi_type< bool >() { return MPI_CXX_BOOL; }
//---------------------------------------------------------------------------
inline void AsterMPI::all_gather( const std::string &in_values, VectorString &out_values,
                                  aster_comm_t *_commCurrent ) {
    // Get number of procs
    const std::size_t comm_size = getMPISize();

    // Get data size on each process
    VectorInt pcounts;
    AsterMPI::all_gather( int( in_values.size() ), pcounts, _commCurrent );

    // Build offsets
    VectorInt offsets( comm_size + 1, 0 );
    for ( std::size_t i = 1; i <= comm_size; ++i )
        offsets[i] = offsets[i - 1] + pcounts[i - 1];

    // Gather
    const std::size_t n = std::accumulate( pcounts.begin(), pcounts.end(), 0 );
    std::vector< char > _out( n );
    aster_mpi_allgatherv( const_cast< char * >( in_values.data() ), in_values.size(), MPI_CHAR,
                          _out.data(), pcounts.data(), offsets.data(), MPI_CHAR, _commCurrent );

    // Rebuild
    out_values.resize( comm_size );
    for ( std::size_t p = 0; p < comm_size; ++p ) {
        out_values[p] = std::string( _out.begin() + offsets[p], _out.begin() + offsets[p + 1] );
    }
}
//---------------------------------------------------------------------------
inline void AsterMPI::all_gather( const VectorString &in_values, VectorString &out_values,
                                  aster_comm_t *_commCurrent ) {
    // Get number of procs
    const std::size_t comm_size = getMPISize();

    // Get data size on each process
    VectorInt pcounts_lc, pcounts_gl;
    pcounts_lc.reserve( in_values.size() );
    for ( auto str : in_values ) {
        pcounts_lc.push_back( str.size() );
    }
    AsterMPI::all_gather( pcounts_lc, pcounts_gl, _commCurrent );

    const std::size_t local_size = std::accumulate( pcounts_lc.begin(), pcounts_lc.end(), 0 );
    const std::size_t global_size = std::accumulate( pcounts_gl.begin(), pcounts_gl.end(), 0 );

    // Gather
    std::vector< char > _in, _out;
    _in.reserve( local_size );
    for ( auto &str : in_values ) {
        for ( auto &carac : str )
            _in.push_back( carac );
    }
    AS_ASSERT( _in.size() == local_size );

    AsterMPI::all_gather( _in, _out, _commCurrent );

    // Rebuild
    VectorInt offsets_gl( pcounts_gl.size() + 1, 0 );
    for ( std::size_t i = 1; i <= pcounts_gl.size(); ++i )
        offsets_gl[i] = offsets_gl[i - 1] + pcounts_gl[i - 1];

    out_values.clear();
    out_values.resize( pcounts_gl.size() );
    for ( std::size_t p = 0; p < pcounts_gl.size(); ++p ) {
        out_values[p] =
            std::string( _out.begin() + offsets_gl[p], _out.begin() + offsets_gl[p + 1] );
    }
}
//---------------------------------------------------------------------------
template < int length >
inline void AsterMPI::all_gather( const std::vector< JeveuxString< length > > &in_values,
                                  std::vector< JeveuxString< length > > &out_values,
                                  aster_comm_t *_commCurrent ) {
    // Get number of procs
    const std::size_t comm_size = getMPISize();
    const std::size_t comm_rank = getMPIRank();

    typedef JeveuxString< length > JeveuxChar;
    VectorInt sizes2;
    AsterMPI::all_gather( int( in_values.size() ), sizes2, _commCurrent );

    out_values.clear();
    for ( int rank = 0; rank < comm_size; ++rank ) {
        JeveuxChar *retour = new JeveuxChar[sizes2[rank]];

        if ( rank == comm_rank )
            for ( int position = 0; position < sizes2[rank]; ++position )
                retour[position] = in_values[position];

        aster_mpi_bcast( retour, length * sizes2[rank], MPI_CHAR, rank, _commCurrent );
        for ( int position = 0; position < sizes2[rank]; ++position )
            out_values.push_back( JeveuxChar( retour[position] ) );
        delete[] retour;
    }
}
//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::all_gather( const std::vector< T > &in_values,
                           std::vector< std::vector< T > > &out_values,
                           aster_comm_t *_commCurrent ) {

    // Get number of procs
    const std::size_t comm_size = getMPISize();

    // Get data size on each process
    VectorInt pcounts;
    const int local_size = in_values.size();
    AsterMPI::all_gather( local_size, pcounts, _commCurrent );
    assert( pcounts.size() == comm_size );

    // Build offsets
    VectorInt offsets( comm_size + 1, 0 );
    for ( std::size_t i = 1; i <= comm_size; ++i )
        offsets[i] = offsets[i - 1] + pcounts[i - 1];

    // Gather data
    const std::size_t n = std::accumulate( pcounts.begin(), pcounts.end(), 0 );
    std::vector< T > recvbuf( n );
    aster_mpi_allgatherv( const_cast< T * >( in_values.data() ), in_values.size(), mpi_type< T >(),
                          recvbuf.data(), pcounts.data(), offsets.data(), mpi_type< T >(),
                          _commCurrent );

    // Repack data
    out_values.resize( comm_size );
    for ( std::size_t p = 0; p < comm_size; ++p ) {
        out_values[p].resize( pcounts[p] );
        for ( int i = 0; i < pcounts[p]; ++i )
            out_values[p][i] = recvbuf[offsets[p] + i];
    }
}
//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::all_gather( const std::vector< T > &in_values, std::vector< T > &out_values,
                           aster_comm_t *_commCurrent ) {

    // Get number of procs
    const std::size_t comm_size = getMPISize();

    // Get data size on each process
    VectorInt pcounts;
    const int local_size = in_values.size();
    AsterMPI::all_gather( local_size, pcounts, _commCurrent );
    assert( pcounts.size() == comm_size );

    // Build offsets
    VectorInt offsets( comm_size + 1, 0 );
    for ( std::size_t i = 1; i <= comm_size; ++i )
        offsets[i] = offsets[i - 1] + pcounts[i - 1];

    // Gather data
    const std::size_t n = std::accumulate( pcounts.begin(), pcounts.end(), 0 );
    out_values.resize( n );
    aster_mpi_allgatherv( const_cast< T * >( in_values.data() ), in_values.size(), mpi_type< T >(),
                          out_values.data(), pcounts.data(), offsets.data(), mpi_type< T >(),
                          _commCurrent );
}
//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::all_gather( const std::vector< T > &in_values, JeveuxVector< T > &out_values,
                           aster_comm_t *_commCurrent ) {

    // Get number of procs
    const std::size_t comm_size = getMPISize();

    // Get data size on each process
    VectorInt pcounts;
    const int local_size = in_values.size();
    AsterMPI::all_gather( local_size, pcounts, _commCurrent );
    assert( pcounts.size() == comm_size );

    // Build offsets
    VectorInt offsets( comm_size + 1, 0 );
    for ( std::size_t i = 1; i <= comm_size; ++i )
        offsets[i] = offsets[i - 1] + pcounts[i - 1];

    // Gather data
    const std::size_t n = std::accumulate( pcounts.begin(), pcounts.end(), 0 );

    if ( out_values->size() < n ) {
        if ( out_values->isAllocated() )
            out_values->deallocate();
        out_values->allocate( n );
    }
    out_values->updateValuePointer();

    aster_mpi_allgatherv( const_cast< T * >( in_values.data() ), in_values.size(), mpi_type< T >(),
                          out_values->getDataPtr(), pcounts.data(), offsets.data(), mpi_type< T >(),
                          _commCurrent );
}
//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::all_gather( const T in_value, std::vector< T > &out_values,
                           aster_comm_t *_commCurrent ) {
    out_values.resize( getMPISize() );
    aster_mpi_allgather( const_cast< T * >( &in_value ), 1, mpi_type< T >(), out_values.data(), 1,
                         mpi_type< T >(), _commCurrent );
}
//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::all_reduce( const T in_value, T &out_value, MPI_Op op, aster_comm_t *_commCurrent ) {
    aster_mpi_allreduce( const_cast< T * >( &in_value ), static_cast< T * >( &out_value ), 1,
                         mpi_type< T >(), op, _commCurrent );
}
//---------------------------------------------------------------------------
template < typename T > void AsterMPI::bcast( T &value, int root, aster_comm_t *_commCurrent ) {
    aster_mpi_bcast( const_cast< T * >( &value ), 1, mpi_type< T >(), root, _commCurrent );
}
//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::bcast( std::vector< T > &value, int root, aster_comm_t *_commCurrent ) {
    aster_mpi_bcast( const_cast< T * >( &value.data() ), value.size(), mpi_type< T >(), root,
                     _commCurrent );
}
//---------------------------------------------------------------------------

#endif /* ASTER_HAVE_MPI */

#endif /* ASTERMPI_H_ */
