#ifndef MEDVECTOR_H
#define MEDVECTOR_H

/**
 * @file MedVector.h
 * @brief Fichier entete de la classe MedVector
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

#include "aster_pybind.h"

#include "MemoryManager/NumpyAccess.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

// aslint: disable=C3010
// aslint: disable=C3008
// aslint: disable=C3012

/**
 * @class MedVector
 * @brief Med vector: 2D vector of size nbelem (nodes or cells) * nbcmp (component or
 * component*nbpg)
 * @author Nicolas Sellenet
 */
template < typename TypeName = double >
class MedVector {
  public:
    typedef TypeName value_type;
    class ElementValue {
        std::vector< TypeName > &_vector;
        int _index = -1;
        int _nbCmp = 0;

        ElementValue( std::vector< TypeName > &vector ) : _vector( vector ) {};

      public:
        int getComponentNumber() const { return _nbCmp; };

        TypeName &operator[]( const int &cmp ) { return _vector[_index + cmp]; };
        const TypeName &operator[]( const int &cmp ) const { return _vector[_index + cmp]; };

        friend class MedVector;
    };

  private:
    /** @brief value vector */
    std::vector< TypeName > _vector;
    /** @brief cumulated size vector (used to explore _vector element by element) */
    std::vector< int > _cumSize;
    /** @brief component number of each element */
    std::vector< int > _cmps;
    /** @brief element number */
    int _size = 0;
    /** @brief component number */
    int _cmpNb = 0;
    /** @brief iterator on current element */
    ElementValue _curVal;
    /** @brief component name */
    std::vector< std::string > _cmpName;

  public:
    /**
     * @brief Constructor
     */
    MedVector( int nbElement )
        : _cumSize( std::vector< int >( nbElement + 1, 0 ) ),
          _cmps( std::vector< int >( nbElement, 0 ) ),
          _size( nbElement ),
          _curVal( ElementValue( _vector ) ) {};

    MedVector() : _curVal( ElementValue( _vector ) ) {};

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    MedVector( const py::tuple &tup ) : MedVector() {
        _vector = tup[0].cast< std::vector< TypeName > >();
        _cumSize = tup[1].cast< std::vector< int > >();
        _cmps = tup[2].cast< std::vector< int > >();
        _cmpNb = tup[3].cast< int >();
        _cmpName = tup[4].cast< std::vector< std::string > >();
    };
    py::tuple _getState() const {
        return py::make_tuple( _vector, _cumSize, _cmps, _cmpNb, _cmpName );
    };

    /** @brief end vector definition -> vector allocation */
    void endDefinition() {
        _cumSize[0] = 0;
        for ( int i = 0; i < _size; ++i ) {
            _cumSize[i + 1] = _cumSize[i] + _cmps[i];
        }
        _vector = std::vector< TypeName >( _cumSize[_size], 0 );
    };

    /** @brief get component name */
    const std::vector< std::string > &getComponentName() const { return _cmpName; };

    /** @brief get component number */
    int getComponentNumber() const { return _cmpNb; };

    /** @brief get component on element vector */
    std::vector< int > getComponentVector() const { return _cmps; };

    /** @brief get cumulated sizes vector */
    std::vector< int > getCumulatedSizesVector() const { return _cumSize; };

    /** @brief get element component number */
    int getElement( int index ) const { return _cmps[index]; };

    /**
     * @brief get values (WARNING values are owned by MedVector: no copy)
     * @return numpy array
     */
    py::object getValues() const {
        npy_intp dims[1] = { (long int)_vector.size() };

        PyObject *values =
            PyArray_SimpleNewFromData( 1, dims, npy_type< double >::value, (void *)&_vector[0] );
        AS_ASSERT( values != NULL );

        PyArray_CLEARFLAGS( (PyArrayObject *)values, NPY_ARRAY_WRITEABLE );
        PyArray_CLEARFLAGS( (PyArrayObject *)values, NPY_ARRAY_OWNDATA );

        py::object tuple = py::reinterpret_steal< py::object >( values );
        tuple.inc_ref();
        return tuple;
    };

    /** @brief set component name */
    void setComponentName( const std::vector< std::string > &cmpName ) {
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

    /** @brief set component name */
    void setComponentVector( const std::vector< int > &cmpVec ) { _cmps = cmpVec; };

    /** @brief set cumulated sizes vector */
    void setCumulatedSizesVector( const std::vector< int > cumSizes ) { _cumSize = cumSizes; };

    /** @brief set element component number */
    void setElement( int index, int nbCmp ) {
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
    void setElements( int index, int size, int nbCmp ) {
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
    void setSize( int nbElement ) {
        _cumSize = std::vector< int >( nbElement + 1, -1 );
        _cmps = std::vector< int >( nbElement, _cmpNb );
        _size = nbElement;
    };

    /** @brief set total size (usefull for memory preallocation) */
    void setTotalSize( int totalSize ) { _vector = std::vector< TypeName >( totalSize ); };

    /** @brief set total element number */
    void setValues( const std::vector< double > &values ) { _vector = values; };

    /** @brief get element number */
    int size() const { return _size; };

    /** @brief get value vector size */
    int totalSize() const { return _vector.size(); };

    const ElementValue &operator[]( const int &index ) const {
        const_cast< MedVector * >( this )->_curVal._index = _cumSize[index];
        const_cast< MedVector * >( this )->_curVal._nbCmp = _cmps[index];
        return _curVal;
    };

    ElementValue &operator[]( const int &index ) {
        _curVal._index = _cumSize[index];
        _curVal._nbCmp = _cmps[index];
        return _curVal;
    };
};

/**
 * @typedef MedVectorPtr
 * @brief Pointeur intelligent vers un MedVector
 */
typedef std::shared_ptr< MedVector< double > > MedVectorPtr;
typedef std::shared_ptr< MedVector< long int > > MedVectorLongPtr;

#endif /* MEDVECTOR_H */
