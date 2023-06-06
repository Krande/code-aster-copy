#ifndef SIMPLEFIELDONNODES_H_
#define SIMPLEFIELDONNODES_H_

/**
 * @file SimpleFieldOnNodes.h
 * @brief Fichier entete de la classe SimpleFieldOnNodes
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

#include "aster_fort_ds.h"

#include "DataFields/DataField.h"
#include "MemoryManager/JeveuxVector.h"
#include "MemoryManager/NumpyAccess.h"
#include "Meshes/BaseMesh.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Tools.h"

/**
 * @class SimpleFieldOnNodes
 * @brief Cette classe template permet de definir un champ aux noeuds Aster
 * @author Nicolas Sellenet
 */
template < class ValueType >
class SimpleFieldOnNodes : public DataField {
  private:
    /** @brief Vecteur Jeveux '.CNSK' */
    JeveuxVectorChar8 _descriptor;
    /** @brief Vecteur Jeveux '.CNSD' */
    JeveuxVectorLong _size;
    /** @brief Vecteur Jeveux '.CNSC' */
    JeveuxVectorChar8 _component;
    /** @brief Vecteur Jeveux '.CNSV' */
    JeveuxVector< ValueType > _values;
    /** @brief Vecteur Jeveux '.CNSL' */
    JeveuxVectorLogical _allocated;
    /** @brief Nombre de noeuds */
    ASTERINTEGER _nbNodes;
    /** @brief Nombre de composantes */
    ASTERINTEGER _nbComp;
    /** @brief Mesh */
    BaseMeshPtr _mesh;

    std::map< std::string, ASTERINTEGER > _name2Index;

    void _buildComponentsName2Index() {
        if ( _name2Index.empty() ) {

            auto nbCmp = this->getNumberOfComponents();
            for ( ASTERINTEGER i = 0; i < nbCmp; i++ ) {
                _name2Index[this->getComponent( i )] = i;
            }
        }
    }

    /**
     * Functions to check an out-of-range condition
     */
    void _checkNodeOOR( const ASTERINTEGER &ino ) const {
        ASTERINTEGER nbNodes = this->getNumberOfNodes();
        if ( ino < 0 || ino >= nbNodes ) {
            throw std::runtime_error( "Node index " + std::to_string( ino ) +
                                      " is out of range [0, " + std::to_string( nbNodes - 1 ) +
                                      "]" );
        };
    }

    void _checkCmpOOR( const ASTERINTEGER &icmp ) const {
        ASTERINTEGER ncmp = this->getNumberOfComponents();
        if ( icmp < 0 || icmp >= ncmp ) {
            throw std::runtime_error( "Component " + std::to_string( icmp ) +
                                      " is out of range [0, " + std::to_string( ncmp - 1 ) + "]" );
        }
    }

    void _checkSize( const ASTERINTEGER &ino, const ASTERINTEGER &icmp ) const {
        if ( this->getNumberOfNodes() == 0 || this->getNumberOfComponents() == 0 )
            throw std::runtime_error( "First call of updateValuePointers is mandatory" );

        this->_checkNodeOOR( ino );
        this->_checkCmpOOR( icmp );
    }

    void _existsValue() const {
        if ( this->getNumberOfNodes() == 0 || this->getNumberOfComponents() == 0 )
            throw std::runtime_error( "First call of updateValuePointers is mandatory" );

        auto size = _allocated->size();
        bool isAlloc = false;
        for ( ASTERINTEGER i = 0; i < size; i++ ) {
            if ( ( *_allocated )[i] ) {
                isAlloc = true;
                break;
            }
        }

        if ( !isAlloc ) {
            throw std::runtime_error( "No values affected in the field" );
        }
    }

  public:
    /**
     * @typedef SimpleFieldOnNodesPtr
     * @brief Pointeur intelligent vers un SimpleFieldOnNodes
     */
    typedef std::shared_ptr< SimpleFieldOnNodes > SimpleFieldOnNodesPtr;

    /**
     * @brief Constructeur
     * @param name Nom Jeveux du champ aux noeuds
     */
    SimpleFieldOnNodes( const std::string name )
        : DataField( name, "CHAM_NO_S" ),
          _descriptor( JeveuxVectorChar8( getName() + ".CNSK" ) ),
          _size( JeveuxVectorLong( getName() + ".CNSD" ) ),
          _component( JeveuxVectorChar8( getName() + ".CNSC" ) ),
          _values( JeveuxVector< ValueType >( getName() + ".CNSV" ) ),
          _allocated( JeveuxVectorLogical( getName() + ".CNSL" ) ),
          _nbNodes( 0 ),
          _nbComp( 0 ) {};

    /**
     * @brief Constructeur

     */
    SimpleFieldOnNodes() : SimpleFieldOnNodes( DataStructureNaming::getNewName( 19 ) ) {};

    SimpleFieldOnNodes( const BaseMeshPtr mesh ) : SimpleFieldOnNodes() { _mesh = mesh; };

    SimpleFieldOnNodes( const BaseMeshPtr mesh, const std::string quantity,
                        const VectorString &comp, bool zero = false )
        : SimpleFieldOnNodes( mesh ) {
        ASTERINTEGER nbComp = comp.size();
        std::string base = "G";

        char *tabNames = vectorStringAsFStrArray( comp, 8 );

        CALL_CNSCRE_WRAP( _mesh->getName().c_str(), quantity.c_str(), &nbComp, tabNames,
                          base.c_str(), getName().c_str(), (ASTERLOGICAL *)&zero );

        FreeStr( tabNames );

        build();
    }

    BaseMeshPtr getMesh() const { return _mesh; };

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    inline ValueType &operator[]( const ASTERINTEGER &i ) { return _values->operator[]( i ); };

    inline const ValueType &operator[]( const ASTERINTEGER &i ) const {
        return _values->operator[]( i );
    };

    ValueType &operator()( const ASTERINTEGER &ino, const ASTERINTEGER &icmp ) {
#ifdef ASTER_DEBUG_CXX
        _checkSize( ino, icmp );
#endif

        const ASTERINTEGER position = ino * this->getNumberOfComponents() + icmp;

        ( *_allocated )[position] = true;
        return this->operator[]( position );
    };

    const ValueType &operator()( const ASTERINTEGER &ino, const ASTERINTEGER &icmp ) const {
#ifdef ASTER_DEBUG_CXX

        if ( !this->hasValue( ino, icmp ) ) {
            AS_ABORT( "DEBUG: Position (" + std::to_string( ino ) + ", " + std::to_string( icmp ) +
                      ") is valid but not allocated!" )
        };
#endif

        const ASTERINTEGER position = ino * this->getNumberOfComponents() + icmp;

        return this->operator[]( position );
    };

    ValueType &operator()( const ASTERINTEGER &ino, const std::string &cmp ) {
        auto icmp = _name2Index.at( cmp );
        return this->operator()( ino, icmp );
    };

    const ValueType &operator()( const ASTERINTEGER &ino, const std::string &cmp ) const {
        auto icmp = _name2Index.at( cmp );
        return this->operator()( ino, icmp );
    };

    bool hasValue( const ASTERINTEGER &ino, const ASTERINTEGER &icmp ) const {

#ifdef ASTER_DEBUG_CXX
        _checkSize( ino, icmp );
#endif

        const ASTERINTEGER position = ino * this->getNumberOfComponents() + icmp;

        return ( *_allocated )[position];
    };

    bool hasValue( const ASTERINTEGER &ino, const std::string &cmp ) const {
        auto icmp = _name2Index.at( cmp );

        return this->hasValue( ino, icmp );
    };

    /**
     * @brief Get number of components
     */
    inline ASTERINTEGER getNumberOfComponents() const { return _nbComp; }

    /**
     * @brief Get number of nodes
     */
    inline ASTERINTEGER getNumberOfNodes() const { return _nbNodes; }

    /**
     * @brief Return a pointer to the vector of data
     */
    const ValueType *getDataPtr() const { return _values->getDataPtr(); }

    /**
     * @brief Get values with mask
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

            Py_DECREF( values );
            Py_DECREF( mask );

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
     * @brief Get the name of the i-th component
     */
    std::string getComponent( const ASTERINTEGER &i ) const {

        if ( i < 0 || i >= _nbComp ) {
            throw std::runtime_error( "Out of range" );
        };

        std::string name = strip( ( *_component )[i].toString() );
        return name;
    };

    /**
     * @Brief Get the names of all the components
     */
    VectorString getComponents() const {

        ASTERINTEGER size = this->getNumberOfComponents();
        VectorString names;
        names.reserve( size );
        for ( ASTERINTEGER i = 0; i < size; i++ ) {
            names.push_back( this->getComponent( i ) );
        }
        return names;
    }

    /**
     * @brief Maps between name of components and the nimber
     */
    std::map< std::string, ASTERINTEGER > getComponentsName2Index() const {

        if ( _name2Index.empty() ) {
            raiseAsterError( "getComponentsName2Index: build field before" );
        }

        return _name2Index;
    };

    /**
     * @brief Get physical quantity
     */
    std::string getPhysicalQuantity() const { return strip( ( *_descriptor )[1].toString() ); }

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    void updateValuePointers() const {
        _descriptor->updateValuePointer();
        _size->updateValuePointer();
        _component->updateValuePointer();
        _values->updateValuePointer();
        _allocated->updateValuePointer();
    };

    bool build() {
        updateValuePointers();

        _nbNodes = ( *_size )[0];
        _nbComp = ( *_size )[1];

        _buildComponentsName2Index();

        AS_ASSERT( _values->size() == _nbNodes * _nbComp );
        AS_ASSERT( _values->size() > 0 );

        return true;
    }

    SimpleFieldOnNodesPtr restrict( const VectorString &cmps = {},
                                    const VectorString &groupsOfNodes = {} ) const {

        this->updateValuePointers();

        VectorString list_cmp;
        auto list_cmp_in = this->getComponents();
        if ( cmps.empty() ) {
            list_cmp = list_cmp_in;
        } else {
            auto set_cmps = toSet( cmps );

            for ( auto &cmp : list_cmp_in ) {
                if ( set_cmps.count( cmp ) > 0 ) {
                    list_cmp.push_back( cmp );
                }
            }
        }

        if ( list_cmp.empty() ) {
            raiseAsterError( "Restriction on list of components is empty" );
        }

        auto ret = std::make_shared< SimpleFieldOnNodes< ValueType > >(
            this->getMesh(), this->getPhysicalQuantity(), list_cmp );

        VectorLong nodes = _mesh->getNodes( groupsOfNodes );

        if ( nodes.empty() ) {
            raiseAsterError( "Restriction on list of nodes is empty" );
        }

        for ( auto &cmp : list_cmp ) {
            auto icmp_in = ( *this )._name2Index.at( cmp );
            auto icmp = ( *ret )._name2Index.at( cmp );
            for ( auto &node : nodes ) {
                if ( this->hasValue( node, icmp_in ) ) {
                    ( *ret )( node, icmp ) = ( *this )( node, icmp_in );
                }
            }
        }

        return ret;
    };
};

/** @typedef SimpleFieldOnNodesReal Class d'une champ simple de doubles */
typedef SimpleFieldOnNodes< ASTERDOUBLE > SimpleFieldOnNodesReal;

/**
 * @typedef SimpleFieldOnNodesPtrReal
 * @brief Definition d'un champ simple de doubles
 */
typedef std::shared_ptr< SimpleFieldOnNodesReal > SimpleFieldOnNodesRealPtr;

/** @typedef SimpleFieldOnNodesLong Class d'un champ simple de long */
typedef SimpleFieldOnNodes< long > SimpleFieldOnNodesLong;

/**
 * @typedef SimpleFieldOnNodesPtrLong
 * @brief Definition d'un champ simple de long
 */
typedef std::shared_ptr< SimpleFieldOnNodesLong > SimpleFieldOnNodesLongPtr;

/** @typedef SimpleFieldOnNodesComplex
    @brief Class d'un champ simple de complexes */
typedef SimpleFieldOnNodes< ASTERCOMPLEX > SimpleFieldOnNodesComplex;

/**
 * @typedef SimpleFieldOnNodesComplexPtr
 * @brief Definition d'un champ simple aux noeuds de complexes
 */
typedef std::shared_ptr< SimpleFieldOnNodesComplex > SimpleFieldOnNodesComplexPtr;
#endif /* SIMPLEFIELDONNODES_H_ */
