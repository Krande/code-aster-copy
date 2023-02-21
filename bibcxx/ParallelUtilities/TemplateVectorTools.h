#ifndef TEMPLATEVECTORTOOLS_H_
#define TEMPLATEVECTORTOOLS_H_

/**
 * @file ObjectBalancer.h
 * @brief Header of tools to manipulate some vector in templates
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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
#include "aster_mpi.h"
#include "astercxx.h"

#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxString.h"
#include "MemoryManager/JeveuxVector.h"

template < typename T >
int getSize( const std::vector< T > &toCopy ) {
    return toCopy.size();
};

template < typename T >
int getSize( const std::vector< std::vector< T > > &toCopy ) {
    return toCopy.size();
};

template < typename T >
int getSize( const JeveuxVector< T > &toCopy ) {
    return toCopy->size();
};

template < typename T >
int getSize( const JeveuxCollection< T > &toCopy ) {
    return toCopy->size();
};

template < typename T >
int getSize( const JeveuxCollectionObject< T > &toCopy ) {
    return toCopy->size();
};

template < typename T >
int getSize( const std::vector< T > *in ) {
    return 1;
};

template < typename T >
int getTotalSize( const std::vector< T > &toCopy ) {
    return 0;
};

template < typename T >
int getTotalSize( const std::vector< std::vector< T > > &toCopy ) {
    return 0;
};

template < typename T >
int getTotalSize( const JeveuxVector< T > &toCopy ) {
    return 0;
};

template < typename T >
int getTotalSize( const JeveuxCollectionClass< T > &toCopy ) {
    return toCopy.totalSize();
};

template < typename T >
const T *getAddress( const std::vector< T > &toCopy ) {
    return &toCopy[0];
};

template < typename T >
T *getAddress( std::vector< T > &toCopy ) {
    return &toCopy[0];
};

template < typename T >
T *getAddress( const JeveuxVector< T > &toCopy ) {
    return &( *toCopy )[0];
};

template < typename T >
T *getAddress( const JeveuxCollection< T > &toCopy ) {
    return &( *toCopy )[0][0];
};

template < typename T >
void resize( std::vector< T > &toCopy, const int &size ) {
    toCopy.resize( size );
};

template < typename T >
void resize( JeveuxVector< T > &toCopy, const int &size ) {
    toCopy->resize( size );
};

template < typename T >
void resize( JeveuxCollection< T > &toCopy, const int &size ) {
    toCopy->resize( size );
};

template < typename T >
struct ValueType;

template < typename T >
struct ValueType< std::vector< T > > {
    typedef T value_type;
};

template < typename T >
struct ValueType< std::vector< std::vector< T > > > {
    typedef T value_type;
};

template < typename T >
struct ValueType< JeveuxVector< T > > {
    typedef T value_type;
};

// template < typename T >
// typename std::enable_if< std::is_integral< T >::value, T >::type
// T& getOccurence( T& )

template < typename T >
const std::vector< T > *getOccurence( const std::vector< std::vector< T > > &in, const int &iPos ) {
    return &in[iPos];
};

template < typename T >
std::vector< T > *getOccurence( std::vector< std::vector< T > > &in, const int &iPos ) {
    return &in[iPos];
};

template < typename T >
JeveuxCollectionObject< T > getOccurence( const JeveuxCollection< T > &in, const int &iPos ) {
    return ( *in )[iPos];
};

template < typename T >
struct StartPosition;

template < typename T >
struct StartPosition< std::vector< std::vector< T > > > {
    static constexpr int value = 0;
};

template < typename T >
struct StartPosition< JeveuxCollectionClass< T > > {
    static constexpr int value = 1;
};

template < typename T >
void allocate( std::vector< std::vector< T > > &in, const int &size1, const int &size2 ) {
    in = std::vector< std::vector< T > >( size1 );
};

template < typename T >
void allocate( std::vector< T > &in, const int &size1, const int &size2 ) {
    in = std::vector< T >( size1 );
};

template < typename T >
void allocate( JeveuxCollectionClass< T > &in, const int &size1, const int &size2 ) {
    in.allocateContiguousNumbered( size1, size2 );
};

template < typename T >
void update( const std::vector< T > &in ){};

template < typename T >
void update( const std::vector< T > *in ){};

template < typename T >
void update( JeveuxCollectionObject< T > in ) {
    in->updateValuePointer();
};

template < typename T >
void allocateOccurence( std::vector< std::vector< T > > &in, const int &pos, const int &size ) {
    in[pos] = std::vector< T >( size );
};

template < typename T >
void allocateOccurence( JeveuxCollectionClass< T > &in, const int &pos, const int &size ) {
    in.allocateObject( pos, size );
};

template < typename T >
const T &getValue( const std::vector< T > &in, const int &pos ) {
    return in[pos];
};

template < typename T >
const T getValue( const std::vector< T > *in, const int &pos ) {
    return ( *in )[pos];
};

template < typename T >
const T &getValue( const JeveuxCollectionObject< T > &in, const int &pos ) {
    return ( *in )[pos];
};

template < typename T >
void setValue( std::vector< T > *in, const int &pos, const T &val ) {
    ( *in )[pos] = val;
};

template < typename T >
void setValue( JeveuxCollectionObject< T > &in, const int &pos, const T &val ) {
    ( *in )[pos] = val;
};

#endif /* TEMPLATEVECTORTOOLS_H_ */
