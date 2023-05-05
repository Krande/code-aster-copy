#ifndef SIMPLEFIELDONCELLS_H_
#define SIMPLEFIELDONCELLS_H_

/**
 * @file SimpleFieldOnCells.h
 * @brief Fichier entete de la classe SimpleFieldOnCells
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

#include "astercxx.h"

#include "DataFields/DataField.h"
#include "MemoryManager/JeveuxVector.h"
#include "MemoryManager/NumpyAccess.h"
#include "Utilities/Tools.h"

/**
 * @class SimpleFieldOnCells
 * @brief Cette classe template permet de definir un champ aux éléments Aster
 * @author Nicolas Sellenet
 */
template < class ValueType >
class SimpleFieldOnCells : public DataField {
  private:
    /** @brief Vecteur Jeveux '.CESK' */
    JeveuxVectorChar8 _descriptor;
    /** @brief Vecteur Jeveux '.CESD' */
    JeveuxVectorLong _size;
    /** @brief Vecteur Jeveux '.CESC' */
    JeveuxVectorChar8 _component;
    /** @brief Vecteur Jeveux '.CESV' */
    JeveuxVector< ValueType > _values;
    /** @brief Vecteur Jeveux '.CESL' */
    JeveuxVectorLogical _allocated;
    /** @brief Nombre de éléments */
    ASTERINTEGER _nbCells;
    /** @brief Nombre de composantes */
    ASTERINTEGER _nbComp;
    /** @brief Number of points */
    ASTERINTEGER _nbPt;
    /** @brief Number of subpoints */
    ASTERINTEGER _nbSpt;

    /**
     * Some unsafe functions to access values without checking dimension
     * Their public version add an if statement
     */

    ASTERINTEGER _ptCell( const ASTERINTEGER &ima ) const { return ( *_size )[4 + 4 * ima + 1]; }
    ASTERINTEGER _sptCell( const ASTERINTEGER &ima ) const { return ( *_size )[4 + 4 * ima + 2]; }
    ASTERINTEGER _cmpsSptCell( const ASTERINTEGER &ima ) const {
        return ( *_size )[4 + 4 * ima + 3];
    }
    ASTERINTEGER _shiftCell( const ASTERINTEGER &ima ) const { return ( *_size )[4 + 4 * ima + 4]; }

    std::string _nameCmp( const ASTERINTEGER &icmp ) const {
        return trim( ( *_component )[icmp].toString() );
    }

    /**
     * Calculate the position of value in CESV array
     */
    ASTERINTEGER _positionInArray( const ASTERINTEGER &icmp, const ASTERINTEGER &ima,
                                   const ASTERINTEGER &ipt, const ASTERINTEGER &ispt ) const {

        ASTERINTEGER npt = this->_ptCell( ima );
        ASTERINTEGER nspt = this->_sptCell( ima );
        ASTERINTEGER ncmp = this->_cmpsSptCell( ima );
        ASTERINTEGER decal = this->_shiftCell( ima );
        return decal + ipt * nspt * ncmp + ispt * ncmp + icmp;
    };

    /**
     * Calculate the size of CESV array
     */
    ASTERINTEGER _nbValArray() const {
        ASTERINTEGER nbVal = 0;
        ASTERINTEGER ncmp_max = this->getNumberOfComponents();
        for ( ASTERINTEGER ima = 0; ima < this->getNumberOfCells(); ima++ ) {
            ASTERINTEGER npt = this->_ptCell( ima );
            ASTERINTEGER nspt = this->_sptCell( ima );
            ASTERINTEGER ncmp = this->_cmpsSptCell( ima );
            if ( ncmp > 0 )
                AS_ASSERT( ncmp == ncmp_max );
            nbVal = nbVal + npt * nspt * ncmp;
        }
        AS_ASSERT( nbVal > 0 );
        return nbVal;
    }

    /**
     * Functions to check an out-of-range condition
     */
    void _checkCellOOR( const ASTERINTEGER &ima ) const {
        ASTERINTEGER nbCells = this->getNumberOfCells();
        if ( ima < 0 || ima >= nbCells ) {
            throw std::runtime_error( "Cell index '" + std::to_string( ima ) +
                                      "' is out of range" );
        };
    }

    void _checkPtOOR( const ASTERINTEGER &ima, const ASTERINTEGER &ipt ) const {
        ASTERINTEGER npt = this->_ptCell( ima );
        if ( ipt < 0 || ipt >= npt ) {
            throw std::runtime_error( "Point     '" + std::to_string( ipt ) +
                                      "' is out of range for cell '" + std::to_string( ima ) +
                                      "'" );
        }
    }

    void _checkSptOOR( const ASTERINTEGER &ima, const ASTERINTEGER &ispt ) const {
        ASTERINTEGER nspt = this->_sptCell( ima );
        if ( ispt < 0 || ispt >= nspt ) {
            throw std::runtime_error( "SubPoint  '" + std::to_string( ispt ) +
                                      "' is out of range for cell '" + std::to_string( ima ) +
                                      "'" );
        }
    }

    void _checkCmpAtCellOOR( const ASTERINTEGER &ima, const ASTERINTEGER &icmp ) const {
        ASTERINTEGER ncmp = this->_cmpsSptCell( ima );
        if ( icmp < 0 || icmp >= ncmp ) {
            throw std::runtime_error( "Component '" + std::to_string( icmp ) +
                                      "' is out of range for cell '" + std::to_string( ima ) +
                                      "'" );
        }
    }

    void _checkCmpOOR( const ASTERINTEGER &icmp ) const {
        ASTERINTEGER ncmp = this->getNumberOfComponents();
        if ( icmp < 0 || icmp >= ncmp ) {
            throw std::runtime_error( "Component '" + std::to_string( icmp ) +
                                      "' is out of range" );
        }
    }

  public:
    /**
     * @typedef SimpleFieldOnCellsPtr
     * @brief Pointeur intelligent vers un SimpleFieldOnCells
     */
    typedef std::shared_ptr< SimpleFieldOnCells > SimpleFieldOnCellsPtr;

    /**
     * @brief Constructeur
     * @param name Nom Jeveux du champ aux éléments
     */
    SimpleFieldOnCells( const std::string name )
        : DataField( name, "CHAM_ELEM_S" ),
          _descriptor( JeveuxVectorChar8( getName() + ".CESK" ) ),
          _size( JeveuxVectorLong( getName() + ".CESD" ) ),
          _component( JeveuxVectorChar8( getName() + ".CESC" ) ),
          _values( JeveuxVector< ValueType >( getName() + ".CESV" ) ),
          _allocated( JeveuxVectorLogical( getName() + ".CESL" ) ),
          _nbCells( 0 ),
          _nbComp( 0 ),
          _nbPt( 0 ),
          _nbSpt( 0 ) {};

    /**
     * @brief Constructeur

     */
    SimpleFieldOnCells() : SimpleFieldOnCells( DataStructureNaming::getNewName( 19 ) ) {};

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    ValueType &operator[]( ASTERINTEGER i ) { return _values->operator[]( i ); };

    /**
     * @brief Access to the (icmp) component of the (ima) cell
              at the (ipt) point, at the (ispt) sub-point.
    */
    ValueType const &getValue( const ASTERINTEGER &ima, const ASTERINTEGER &icmp,
                               const ASTERINTEGER &ipt, const ASTERINTEGER &ispt ) const {

#ifdef ASTER_DEBUG_CXX
        if ( this->getNumberOfCells() == 0 || this->getNumberOfComponents() == 0 )
            throw std::runtime_error( "First call of updateValuePointers is mandatory" );
#endif

        ASTERINTEGER position = this->_positionInArray( icmp, ima, ipt, ispt );

#ifdef ASTER_DEBUG_CXX
        bool allocated = ( *_allocated )[position];
        if ( !allocated ) {
            std::cout << "DEBUG: Position (" + std::to_string( icmp ) + ", " +
                             std::to_string( ima ) + ", " + std::to_string( ipt ) + ", " +
                             std::to_string( ispt ) + ") is valid but not allocated!"
                      << std::endl;
        };
#endif

        return ( *_values )[position];
    }

    /**
     * @brief tell if value exists for (icmp) component of the (ima) cell
              at the (ipt) point, at the (ispt) sub-point.
    */
    bool hasValue( const ASTERINTEGER &ima, const ASTERINTEGER &icmp, const ASTERINTEGER &ipt,
                   const ASTERINTEGER &ispt ) const {

        ASTERINTEGER position = this->_positionInArray( icmp, ima, ipt, ispt );

#ifdef ASTER_DEBUG_CXX
        if ( this->getNumberOfCells() == 0 || this->getNumberOfComponents() == 0 )
            throw std::runtime_error( "First call of updateValuePointers is mandatory" );

        this->_checkCellOOR( ima );
        this->_checkPtOOR( ima, ipt );
        this->_checkSptOOR( ima, ispt );
        this->_checkCmpAtCellOOR( ima, icmp );
#endif

        return ( *_allocated )[position];
    }

    /**
     * @brief Get number of points of the i-th cell
     */
    ASTERINTEGER getNumberOfPointsOfCell( const ASTERINTEGER &ima ) const {
        this->_checkCellOOR( ima );
        return this->_ptCell( ima );
    }

    /**
     * @brief Get number of sub-points of the i-th cell
     */
    ASTERINTEGER getNumberOfSubPointsOfCell( const ASTERINTEGER &ima ) const {
        this->_checkCellOOR( ima );
        return this->_sptCell( ima );
    }

    /**
     * @brief Get number of components of the i-th cell
     */
    ASTERINTEGER getNumberOfComponentsForSubpointsOfCell( const ASTERINTEGER &ima ) const {
        this->_checkCellOOR( ima );
        return this->_cmpsSptCell( ima );
    }

    /**
     * @brief Get number of components
     */
    ASTERINTEGER getNumberOfComponents() const { return _nbComp; }

    /**
     * @brief Get number of nodes
     */
    ASTERINTEGER getNumberOfCells() const { return _nbCells; }

    /**
     * @brief Get number of points
     */
    ASTERINTEGER getMaxNumberOfPoints() const { return _nbPt; }

    /**
     * @brief Get number of sub-points
     */
    ASTERINTEGER getMaxNumberOfSubPoints() const { return _nbSpt; }

    /**
     * @brief Get the name of the i-th component
     */
    std::string getComponent( const ASTERINTEGER &icmp ) const {

        if ( icmp < 0 || icmp >= this->getNumberOfComponents() ) {
            throw std::runtime_error( "Component '" + std::to_string( icmp ) +
                                      "' is out of range" );
        };
        return this->_nameCmp( icmp );
    };

    /**
     * @Brief Get the names of all the components
     */
    VectorString getComponents() const {

        ASTERINTEGER size = this->getNumberOfComponents();
        VectorString names;
        names.reserve( size );
        for ( ASTERINTEGER icmp = 0; icmp < size; icmp++ ) {
            names.push_back( this->_nameCmp( icmp ) );
        }
        return names;
    }

    /**
     * @brief Get physical quantity
     */
    std::string getPhysicalQuantity() const { return trim( ( *_descriptor )[1].toString() ); }

    /**
     * @brief Get field location
     */
    std::string getFieldLocation() const { return trim( ( *_descriptor )[2].toString() ); }

    /**
     * @brief Get cells holding components
     */
    VectorLong getCellsWithComponents() const {
        VectorLong values;
        for ( ASTERINTEGER ima = 0; ima < this->getNumberOfCells(); ima++ ) {
            if ( this->_cmpsSptCell( ima ) > 0 )
                values.push_back( ima );
        }
        return values;
    }

    /**
     * @brief Get values on cells holding components, with mask
     */
    py::object getValues( bool copy = false ) {

        PyObject *resu_tuple = PyTuple_New( 2 );

        npy_intp dims[2] = {_values->size() / this->getNumberOfComponents(),
                            this->getNumberOfComponents()};

        PyObject *values = PyArray_SimpleNewFromData( 2, dims, npy_type< ValueType >::value,
                                                      _values->getDataPtr() );
        PyObject *mask = PyArray_SimpleNewFromData( 2, dims, NPY_BOOL, _allocated->getDataPtr() );
        AS_ASSERT( values != NULL );
        AS_ASSERT( mask != NULL );

        if ( copy ) {
            PyObject *values_copy =
                PyArray_NewLikeArray( (PyArrayObject *)values, NPY_ANYORDER, NULL, 0 );
            PyArray_CopyInto( (PyArrayObject *)values_copy, (PyArrayObject *)values );
            AS_ASSERT( values_copy != NULL );

            PyObject *mask_copy =
                PyArray_NewLikeArray( (PyArrayObject *)mask, NPY_ANYORDER, NULL, 0 );
            PyArray_CopyInto( (PyArrayObject *)mask_copy, (PyArrayObject *)mask );
            AS_ASSERT( mask_copy != NULL );

            PyArray_ENABLEFLAGS( (PyArrayObject *)values_copy, NPY_ARRAY_OWNDATA );
            PyArray_ENABLEFLAGS( (PyArrayObject *)mask_copy, NPY_ARRAY_OWNDATA );

            Py_XDECREF( values );
            Py_XDECREF( mask );

            PyTuple_SetItem( resu_tuple, 0, values_copy );
            PyTuple_SetItem( resu_tuple, 1, mask_copy );

        } else {
            PyArray_CLEARFLAGS( (PyArrayObject *)values, NPY_ARRAY_WRITEABLE );
            PyArray_CLEARFLAGS( (PyArrayObject *)mask, NPY_ARRAY_WRITEABLE );
            PyArray_CLEARFLAGS( (PyArrayObject *)values, NPY_ARRAY_OWNDATA );
            PyArray_CLEARFLAGS( (PyArrayObject *)mask, NPY_ARRAY_OWNDATA );
            PyTuple_SetItem( resu_tuple, 0, values );
            PyTuple_SetItem( resu_tuple, 1, mask );
        }

        py::object tuple = py::reinterpret_steal< py::object >( resu_tuple );
        tuple.inc_ref();
        return tuple;
    }

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    void updateValuePointers() {
        _descriptor->updateValuePointer();
        _size->updateValuePointer();
        _component->updateValuePointer();
        _values->updateValuePointer();
        _allocated->updateValuePointer();
        _nbCells = ( *_size )[0];
        _nbComp = ( *_size )[1];
        _nbPt = ( *_size )[2];
        _nbSpt = ( *_size )[3];

        AS_ASSERT( _values->size() == this->_nbValArray() );
    };
};

using SimpleFieldOnCellsReal = SimpleFieldOnCells< ASTERDOUBLE >;
using SimpleFieldOnCellsRealPtr = std::shared_ptr< SimpleFieldOnCellsReal >;
using SimpleFieldOnCellsLong = SimpleFieldOnCells< ASTERINTEGER >;
using SimpleFieldOnCellsLongPtr = std::shared_ptr< SimpleFieldOnCellsLong >;

#endif /* SIMPLEFIELDONCELLS_H_ */
