#ifndef ARRAYWRAPPER_H_
#define ARRAYWRAPPER_H_

/**
 * @file ArrayWrapper.h
 * @brief Fichier entete de la classe ArrayWrapper
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "aster_pybind.h"

// #include "IOManager/ArrayWrapper< ClassWrap >.h"
#include "MemoryManager/JeveuxVector.h"
#include "MemoryManager/NumpyAccess.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

template < typename T1 >
class VectorInterface;

template < class T >
class VectorInterface< std::vector< T > > {
    /** @brief Smart pointer on object */
    std::vector< T > &_object;

  public:
    typedef T value_type;

    VectorInterface() = delete;

    VectorInterface( std::vector< T > &obj ) : _object( obj ) {};

    void allocate( const int &nbElem ) { _object = std::vector< T >( nbElem, -1 ); };

    void deallocate() { _object = std::vector< T >(); };

    value_type *getPointer() const {
        const T &ref = _object[0];
        T *pointer = const_cast< T * >( &ref );
        return pointer;
    };

    ASTERINTEGER totalSize() const { return _object.size(); };
};

template < class T >
class VectorInterface< JeveuxVector< T > > {
    /** @brief Smart pointer on object */
    JeveuxVector< T > _object;

  public:
    typedef T value_type;

    VectorInterface() = delete;

    VectorInterface( JeveuxVector< T > &obj ) : _object( obj ) {};

    void allocate( const int &nbElem ) { _object->allocate( nbElem ); };

    void deallocate() { _object->deallocate(); };

    value_type *getPointer() const {
        if ( _object->exists() ) {
            _object->updateValuePointer();
            return &( ( *_object )[0] );
        } else {
            return nullptr;
        }
    };

    ASTERINTEGER totalSize() const { return _object->size(); };
};

/**
 * @class ArrayWrapper
 * @brief 2D vector of size nbelem (nodes or cells) * nbcmp (component or
 * component*nbpg)
 * @author Nicolas Sellenet
 */
template < typename ClassWrap >
class ArrayWrapper : public VectorInterface< ClassWrap > {
  public:
    class ElementValue {
        typedef typename VectorInterface< ClassWrap >::value_type value_type;
        value_type *_vector;
        ASTERINTEGER _index = -1;
        ASTERINTEGER _nbCmp = 0;

        ElementValue( value_type *vec ) : _vector( vec ) {};

      protected:
        void setPointer( typename VectorInterface< ClassWrap >::value_type *ptr ) {
            _vector = ptr;
        };

      public:
        ASTERINTEGER getComponentNumber() const { return _nbCmp; };

        value_type &operator[]( const int &cmp ) { return _vector[_index + cmp]; };
        const value_type &operator[]( const int &cmp ) const { return _vector[_index + cmp]; };

        friend class ArrayWrapper< ClassWrap >;
    };

  private:
    /** @brief value vector */
    typename VectorInterface< ClassWrap >::value_type *_pointer;
    /** @brief cumulated size vector (used to explore _vector element by element) */
    VectorLong _cumSize;
    /** @brief component number of each element */
    VectorLong _cmps;
    /** @brief element number */
    ASTERINTEGER _size = 0;
    /** @brief component number */
    ASTERINTEGER _cmpNb = 0;
    /** @brief iterator on current element */
    ElementValue _curVal;
    /** @brief component name */
    VectorString _cmpName;

  public:
    /**
     * @brief Constructors
     */
    ArrayWrapper( ClassWrap &curObj, ASTERINTEGER cmpNb )
        : VectorInterface< ClassWrap >( curObj ),
          _pointer( VectorInterface< ClassWrap >::getPointer() ),
          _cmpNb( cmpNb ),
          _curVal( ElementValue( _pointer ) ) {
        const auto size = VectorInterface< ClassWrap >::totalSize();
        if ( size != 0 ) {
            const int v1 = size / _cmpNb;
            const double v2 = size / _cmpNb;
            if ( v1 != v2 ) {
                throw std::runtime_error( "Inconsistent component number" );
            }
            _size = v1;
            int pos = 0;
            _cmps = VectorLong( _size, _cmpNb );
            _cumSize.reserve( _size + 1 );
            for ( int i = 0; i < _size + 1; ++i ) {
                _cumSize.push_back( pos );
                pos += _cmpNb;
            }
        }
    };

    ArrayWrapper( ClassWrap &curObj, const VectorLong &cmps )
        : VectorInterface< ClassWrap >( curObj ),
          _pointer( VectorInterface< ClassWrap >::getPointer() ),
          _cmps( cmps ),
          _curVal( ElementValue( _pointer ) ) {
        const auto size = VectorInterface< ClassWrap >::totalSize();
        if ( size != 0 ) {
            _size = cmps.size();
            int pos = 0;
            _cmpNb = 0;
            _cumSize.reserve( _size + 1 );
            for ( int i = 0; i < _size; ++i ) {
                const auto &curCmpNb = _cmps[i];
                _cmpNb = std::max( _cmpNb, curCmpNb );
                _cumSize.push_back( pos );
                pos += curCmpNb;
            }
            _cumSize.push_back( pos );
            if ( size != pos ) {
                throw std::runtime_error( "Sizes are not consistents" );
            }
        }
    };

    void deallocate() {
        VectorInterface< ClassWrap >::deallocate();
        _cumSize = VectorLong();
        _cmps = VectorLong();
        _size = 0;
        _cmpNb = 0;
        _cmpName = VectorString();
    };

    /** @brief end vector definition -> vector allocation */
    void endDefinition() {
        _cumSize[0] = 0;
        for ( int i = 0; i < _size; ++i ) {
            _cumSize[i + 1] = _cumSize[i] + _cmps[i];
        }
        VectorInterface< ClassWrap >::allocate( _cumSize[_size] );
        _pointer = VectorInterface< ClassWrap >::getPointer();
        _curVal.setPointer( _pointer );
    };

    /** @brief get component name */
    const VectorString &getComponentName() const { return _cmpName; };

    /** @brief get component number */
    ASTERINTEGER getComponentNumber() const { return _cmpNb; };

    /** @brief get element component number */
    int getElement( int index ) const { return _cmps[index]; };

    /** @brief set component name */
    void setComponentName( const VectorString &cmpName ) {
        if ( _cmpNb == 0 ) {
            _cmpNb = cmpName.size();
        } else {
            if ( cmpName.size() != 0 && cmpName.size() != _cmpNb )
                throw std::runtime_error( "Bad component number" );
        }
        _cmpName = cmpName;
    };

    /** @brief set component number */
    void setComponentNumber( int nbCmp ) { _cmpNb = nbCmp; };

    /** @brief set element component number */
    void setElement( ASTERINTEGER index, ASTERINTEGER nbCmp ) {
        if ( _size == 0 )
            throw std::runtime_error( "Number of elements must be set before set elements" );
        // if ( nbCmp < 1 )
        //     throw std::runtime_error( "Component numbers must be grower than 0" );
        if ( index < 0 || index > _size )
            throw std::runtime_error( "Out of bound" );
        _cmps[index] = nbCmp;
        if ( index == 0 )
            _cumSize[index] = 0;
        if ( index > 0 && _cumSize[index - 1] >= 0 )
            _cumSize[index] = _cumSize[index - 1] + _cmps[index - 1];
    };

    /** @brief set elements component number (slice) */
    void setElements( ASTERINTEGER index, ASTERINTEGER size, ASTERINTEGER nbCmp ) {
        if ( _size == 0 )
            throw std::runtime_error( "Number of elements must be set before set elements" );
        if ( nbCmp < 1 )
            throw std::runtime_error( "Component numbers must be grower than 0" );
        if ( index < 0 || index + size > _size )
            throw std::runtime_error( "Out of bound" );
        for ( int i = 0; i < size; ++i ) {
            _cmps[index + i] = nbCmp;
        }
    };

    /** @brief set total element number */
    void setSize( ASTERINTEGER nbElement ) {
        _cumSize = VectorLong( nbElement + 1, -1 );
        _cmps = VectorLong( nbElement, _cmpNb );
        _size = nbElement;
    };

    /** @brief set total size (usefull for memory preallocation) */
    void setTotalSize( ASTERINTEGER totalSize ) {
        VectorInterface< ClassWrap >::allocate( totalSize );
        _pointer = VectorInterface< ClassWrap >::getPointer();
        _curVal.setPointer( _pointer );
    };

    /** @brief get element number */
    ASTERINTEGER size() const { return _size; };

    /** @brief get value vector size */
    ASTERINTEGER totalSize() const { return VectorInterface< ClassWrap >::totalSize(); };

    const ElementValue &operator[]( const int &index ) const {
        const_cast< ArrayWrapper< ClassWrap > * >( this )->_curVal._index = _cumSize[index];
        const_cast< ArrayWrapper< ClassWrap > * >( this )->_curVal._nbCmp = _cmps[index];
        return _curVal;
    };

    ElementValue &operator[]( const int &index ) {
        _curVal._index = _cumSize[index];
        _curVal._nbCmp = _cmps[index];
        return _curVal;
    };
};

/**
 * @typedef ArrayWrapperPtr
 * @brief Pointeur intelligent vers un ArrayWrapper
 */
typedef std::shared_ptr< ArrayWrapper< JeveuxVector< double > > > ArrayWrapperPtr;
typedef std::shared_ptr< ArrayWrapper< JeveuxVector< bool > > > ArrayWrapperLogicalPtr;
typedef std::shared_ptr< ArrayWrapper< JeveuxVector< ASTERINTEGER > > > ArrayWrapperLongPtr;
typedef std::shared_ptr< ArrayWrapper< VectorReal > > ArrayWrapperVectorReal;
typedef std::shared_ptr< ArrayWrapper< VectorLong > > ArrayWrapperVectorLong;

#endif /* ARRAYWRAPPER_H_ */
