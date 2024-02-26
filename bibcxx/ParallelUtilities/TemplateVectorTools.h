#ifndef TEMPLATEVECTORTOOLS_H_
#define TEMPLATEVECTORTOOLS_H_

/**
 * @file ObjectBalancer.h
 * @brief Header of tools to manipulate some vector in templates
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */
#include "astercxx.h"

#include "aster_mpi.h"

#include "IOManager/MedVector.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxString.h"
#include "MemoryManager/JeveuxVector.h"
#include "ParallelUtilities/ArrayWrapper.h"

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

#ifdef ASTER_HAVE_MED
int getSize( const MedVector< double >::ElementValue &in );
int getSize( const MedVector< long int >::ElementValue &in );
#endif

int getSize( const ArrayWrapper< JeveuxVectorReal >::ElementValue &in );

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

#ifdef ASTER_HAVE_MED
int getTotalSize( const MedVector< double > &toCopy );
int getTotalSize( const MedVector< long int > &toCopy );
#endif

template < typename T >
int getTotalSize( const ArrayWrapper< T > &toCopy ) {
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
    return toCopy->getDataPtr();
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

template < typename T >
struct ValueType< ArrayWrapper< JeveuxVector< T > > > {
    typedef T value_type;
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

#ifdef ASTER_HAVE_MED
template <>
struct StartPosition< MedVector< double > > {
    static constexpr int value = 0;
};
template <>
struct StartPosition< MedVector< long int > > {
    static constexpr int value = 0;
};
#endif
template < typename T >
struct StartPosition< ArrayWrapper< T > > {
    static constexpr int value = 0;
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

#ifdef ASTER_HAVE_MED
void allocate( MedVector< double > &in, const int &size1, const int &size2 );
void allocate( MedVector< long int > &in, const int &size1, const int &size2 );
#endif

template < typename T >
void allocate( ArrayWrapper< T > &in, const int &size1, const int &size2 ) {
    in.setSize( size1 );
    in.setTotalSize( size2 );
};

template < typename T >
void update( const std::vector< T > &in ) {};

template < typename T >
void update( const std::vector< T > *in ) {};

template < typename T >
void update( JeveuxCollectionObject< T > in ) {
    in->updateValuePointer();
};

#ifdef ASTER_HAVE_MED
void update( MedVector< double >::ElementValue in );
void update( MedVector< long int >::ElementValue in );
#endif

void update( typename ArrayWrapper< JeveuxVectorReal >::ElementValue in );

template < typename T >
void allocateOccurence( std::vector< std::vector< T > > &in, const int &pos, const int &size ) {
    in[pos] = std::vector< T >( size );
};

template < typename T >
void allocateOccurence( JeveuxCollectionClass< T > &in, const int &pos, const int &size ) {
    in.allocateObject( pos, size );
};

#ifdef ASTER_HAVE_MED
void allocateOccurence( MedVector< double > &in, const int &pos, const int &size );
void allocateOccurence( MedVector< long int > &in, const int &pos, const int &size );
#endif

template < typename T >
void allocateOccurence( ArrayWrapper< T > &in, const int &pos, const int &size ) {
    in.setElement( pos, size );
};

template < typename T >
struct ObjectTemplateType;

#ifdef ASTER_HAVE_MED
template <>
struct ObjectTemplateType< MedVector< double > > {
    typedef double value_type;
};
template <>
struct ObjectTemplateType< MedVector< long int > > {
    typedef long int value_type;
};
#endif

template <>
struct ObjectTemplateType< JeveuxCollectionClass< ASTERINTEGER > > {
    typedef ASTERINTEGER value_type;
};

template <>
struct ObjectTemplateType< VectorOfVectorsLong > {
    typedef ASTERINTEGER value_type;
};

template < typename T >
struct ObjectTemplateType< ArrayWrapper< JeveuxVector< T > > > {
    typedef T value_type;
};

#endif /* TEMPLATEVECTORTOOLS_H_ */
