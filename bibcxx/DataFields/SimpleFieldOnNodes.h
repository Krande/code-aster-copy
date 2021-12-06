#ifndef SIMPLEFIELDONNODES_H_
#define SIMPLEFIELDONNODES_H_

/**
 * @file SimpleFieldOnNodes.h
 * @brief Fichier entete de la classe SimpleFieldOnNodes
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

#include "MemoryManager/JeveuxVector.h"
#include "DataStructures/DataStructure.h"
#include "Utilities/Tools.h"

/**
 * @class SimpleFieldOnNodes
 * @brief Cette classe template permet de definir un champ aux noeuds Aster
 * @author Nicolas Sellenet
 */
template < class ValueType > class SimpleFieldOnNodes : public DataStructure {
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
    int _nbNodes;
    /** @brief Nombre de composantes */
    int _nbComp;

  public:
    /**
     * @typedef SimpleFieldOnNodesPtr
     * @brief Pointeur intelligent vers un SimpleFieldOnNodes
     */
    typedef boost::shared_ptr< SimpleFieldOnNodes > SimpleFieldOnNodesPtr;

    /**
     * @brief Constructeur
     * @param name Nom Jeveux du champ aux noeuds
     */
    SimpleFieldOnNodes( const std::string name )
        : DataStructure( name, 19, "CHAM_NO_S" ),
          _descriptor( JeveuxVectorChar8( getName() + ".CNSK" ) ),
          _size( JeveuxVectorLong( getName() + ".CNSD" ) ),
          _component( JeveuxVectorChar8( getName() + ".CNSC" ) ),
          _values( JeveuxVector< ValueType >( getName() + ".CNSV" ) ),
          _allocated( JeveuxVectorLogical( getName() + ".CNSL" ) ), _nbNodes( 0 ), _nbComp( 0 ){};

    /**
     * @brief Constructeur

     */
    SimpleFieldOnNodes(  )
        : SimpleFieldOnNodes( DataStructureNaming::getNewName( 19 ) ){};


    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    ValueType &operator[]( int i ) { return _values->operator[]( i ); };

    const ValueType &getValue( int nodeNumber, int compNumber ) const
    {
#ifdef ASTER_DEBUG_CXX
        if ( _nbNodes == 0 || _nbComp == 0 )
            throw std::runtime_error( "First call of updateValuePointers is mandatory" );
#endif

    if ( nodeNumber < 0 || nodeNumber > _nbNodes) {
      throw std::runtime_error("Out of range");
    };
    if ( compNumber < 0 || compNumber > _nbComp) {
      throw std::runtime_error("Out of range");
    };

        const long position = nodeNumber * _nbComp + compNumber;
        return ( *_values )[position];
    };

    /**
     * @brief Get number of components
     */
    int getNumberOfComponents() const { return _nbComp; }

    /**
     * @brief Get number of nodes
     */
    int getNumberOfNodes() const { return _nbNodes; }

    /**
     * @brief Return a pointer to the vector of data
     */
    const ValueType *getDataPtr() const { return _values->getDataPtr(); }

    /**
     * @brief Get values
     */
    std::vector< ValueType > getValues() const {
      const ValueType *ptr = _values->getDataPtr();
      std::vector< ValueType > vect( ptr, ptr + _values->size() );
      return vect;
    }

    /**
     * @brief Get the name of the i-th component
     */
    std::string getNameOfComponent( const int& i ) const {

      if ( i < 0 || i >= _nbComp) {
        throw std::runtime_error( "Out of range");
      };

      std::string name = trim(( *_component )[i].toString());
      return name;
    };


    /**
     * @Brief Get the names of all the components
     */
    VectorString getNameOfComponents() const {

      int size = this->getNumberOfComponents();
      VectorString names;
      names.reserve(size);
      for ( int i = 0 ; i < size; i++ ) {
        names.push_back(this->getNameOfComponent(i));
      }
      return names;
    }

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    bool updateValuePointers() {
        bool retour = _descriptor->updateValuePointer();
        retour = ( retour && _size->updateValuePointer() );
        retour = ( retour && _component->updateValuePointer() );
        retour = ( retour && _values->updateValuePointer() );
        retour = ( retour && _allocated->updateValuePointer() );
        if ( retour ) {
            _nbNodes = ( *_size )[0];
            _nbComp = ( *_size )[1];
            if ( _values->size() != _nbNodes * _nbComp )
                throw std::runtime_error( "Programming error" );
        }
        return retour;
    };
};

/** @typedef SimpleFieldOnNodesReal Class d'une champ simple de doubles */
typedef SimpleFieldOnNodes< ASTERDOUBLE > SimpleFieldOnNodesReal;

/**
 * @typedef SimpleFieldOnNodesPtrReal
 * @brief Definition d'un champ simple de doubles
 */
typedef boost::shared_ptr< SimpleFieldOnNodesReal > SimpleFieldOnNodesRealPtr;

/** @typedef SimpleFieldOnNodesLong Class d'un champ simple de long */
typedef SimpleFieldOnNodes< long > SimpleFieldOnNodesLong;

/**
 * @typedef SimpleFieldOnNodesPtrLong
 * @brief Definition d'un champ simple de long
 */
typedef boost::shared_ptr< SimpleFieldOnNodesLong > SimpleFieldOnNodesLongPtr;

/** @typedef SimpleFieldOnNodesComplex
    @brief Class d'un champ simple de complexes */
typedef SimpleFieldOnNodes< ASTERCOMPLEX > SimpleFieldOnNodesComplex;

/**
 * @typedef SimpleFieldOnNodesComplexPtr
 * @brief Definition d'un champ simple aux noeuds de complexes
 */
typedef boost::shared_ptr< SimpleFieldOnNodesComplex > SimpleFieldOnNodesComplexPtr;
#endif /* SIMPLEFIELDONNODES_H_ */
