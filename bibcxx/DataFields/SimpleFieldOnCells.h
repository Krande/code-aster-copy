#ifndef SIMPLEFIELDONCELLS_H_
#define SIMPLEFIELDONCELLS_H_

/**
 * @file SimpleFieldOnCells.h
 * @brief Fichier entete de la classe SimpleFieldOnCells
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

#include <string>
#include <assert.h>

#include "astercxx.h"

#include "MemoryManager/NumpyAccess.h"
#include "MemoryManager/JeveuxVector.h"
#include "DataStructures/DataStructure.h"
#include "Utilities/Tools.h"

/**
 * @class SimpleFieldOnCells
 * @brief Cette classe template permet de definir un champ aux éléments Aster
 * @author Nicolas Sellenet
 */
template < class ValueType > class SimpleFieldOnCells : public DataStructure {
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
    int _nbCells;
    /** @brief Nombre de composantes */
    int _nbComp;
    /** @brief Number of points */
    int _nbPt;
    /** @brief Number of under points */
    int _nbSpt;

    /**
     * Some unsafe functions to access values without checking dimension
     * Their public version add an if statement
     */

    int _ptCell( const int& ima )      const { return ( *_size )[4 + 4*ima + 1]; }
    int _sptCell( const int& ima )     const { return ( *_size )[4 + 4*ima + 2]; }
    int _cmpsSptCell( const int& ima ) const { return ( *_size )[4 + 4*ima + 3]; }
    int _shiftCell( const int& ima )   const { return ( *_size )[4 + 4*ima + 4]; }

    std::string _nameCmp( const int& icmp) const { return trim(( *_component )[icmp].toString()); }

    /**
     * Calculate the position of value in CESV array
     */
    int _positionInArray( const int& icmp, const int& ima,
                 const int& ipt, const int& ispt) const {

      int npt = this->_ptCell(ima);
      int nspt = this->_sptCell(ima);
      int ncmp = this->_cmpsSptCell(ima);
      int decal = this->_shiftCell(ima);
      return decal + ipt*nspt*ncmp + ispt*ncmp + icmp;
    };

    /**
     * Calculate the size of CESV array
     */
    int _nbValArray( ) const {
      int nbVal = 0;
      for (int ima=0 ; ima<this->getNumberOfCells() ;  ima++){
        int npt = this->_ptCell(ima);
        int nspt = this->_sptCell(ima);
        int ncmp = this->_cmpsSptCell(ima);
        nbVal = nbVal + npt*nspt*ncmp;
      }
      return nbVal;
    }

    /**
     * Functions to check an out-of-range condition
     */
    void _checkCellOOR( const int& ima ) const {
      int nbCells = this->getNumberOfCells();
      if ( ima < 0 || ima >= nbCells) {
        throw std::runtime_error( "Cell index '"+std::to_string(ima)+"' is out of range");
      };
    }

    void _checkPtOOR( const int& ima, const int& ipt ) const {
      int npt = this->_ptCell(ima);
      if ( ipt  < 0 || ipt  >= npt  ) {
        throw std::runtime_error( "Point     '"+std::to_string(ipt)
                                  + "' is out of range for cell '"
                                  +std::to_string(ima)+"'");
      }
    }

    void _checkSptOOR( const int& ima, const int& ispt ) const {
      int nspt = this->_sptCell(ima);
      if ( ispt < 0 || ispt >= nspt ) {
        throw std::runtime_error( "SubPoint  '"+std::to_string(ispt)
                                  +"' is out of range for cell '"
                                  +std::to_string(ima)+"'");
      }
    }

    void _checkCmpAtCellOOR( const int& ima, const int& icmp ) const {
      int ncmp = this->_cmpsSptCell(ima);
      if ( icmp < 0 || icmp >= ncmp ) {
        throw std::runtime_error( "Component '"+std::to_string(icmp)
                                  +"' is out of range for cell '"
                                  +std::to_string(ima)+"'");
      }
    }

    void _checkCmpOOR( const int& icmp ) const {
      int ncmp = this->getNumberOfComponents();
      if ( icmp < 0 || icmp >= ncmp ) {
        throw std::runtime_error( "Component '"+std::to_string(icmp)
                                  +"' is out of range");
      }
    }

  public:
    /**
     * @typedef SimpleFieldOnCellsPtr
     * @brief Pointeur intelligent vers un SimpleFieldOnCells
     */
    typedef boost::shared_ptr< SimpleFieldOnCells > SimpleFieldOnCellsPtr;

    /**
     * @brief Constructeur
     * @param name Nom Jeveux du champ aux éléments
     */
    SimpleFieldOnCells( const std::string name )
        : DataStructure( name, 19, "CHAM_ELEM_S" ),
          _descriptor( JeveuxVectorChar8( getName() + ".CESK" ) ),
          _size( JeveuxVectorLong( getName() + ".CESD" ) ),
          _component( JeveuxVectorChar8( getName() + ".CESC" ) ),
          _values( JeveuxVector< ValueType >( getName() + ".CESV" ) ),
          _allocated( JeveuxVectorLogical( getName() + ".CESL" ) ),
          _nbCells( 0 ), _nbComp( 0 ), _nbPt( 0 ), _nbSpt( 0 ){};

    /**
     * @brief Constructeur

     */
    SimpleFieldOnCells(  )
        : SimpleFieldOnCells( DataStructureNaming::getNewName( 19 ) ){};

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    ValueType &operator[]( int i ) { return _values->operator[]( i ); };

    /**
     * @brief Access to the (icmp) component of the (ima) cell
              at the (ipt) point, at the (ispt) sub-point.
    */
    ValueType const &getValue( const int& ima, const int& icmp,
                               const int& ipt, const int& ispt  ) const {

#ifdef ASTER_DEBUG_CXX
      if ( this->getNumberOfCells() == 0 || this->getNumberOfComponents() == 0)
        throw std::runtime_error( "First call of updateValuePointers is mandatory" );
#endif

      this->_checkCellOOR(ima);
      this->_checkPtOOR(ima, ipt);
      this->_checkSptOOR(ima, ispt);
      this->_checkCmpAtCellOOR(ima, icmp);

      int position = this->_positionInArray(icmp, ima, ipt, ispt);
      bool allocated = ( *_allocated )[position];

#ifdef ASTER_DEBUG_CXX
      if ( !allocated ){
        std::cout <<
          "DEBUG: Position ("+std::to_string(icmp)
          +", "+std::to_string(ima)
          +", "+std::to_string(ipt)
          +", "+std::to_string(ispt)
          +") is valid but not allocated!"
                  << std::endl;
      };
#endif

      return ( *_values )[position];
    }

    /**
     * @brief Get number of points of the i-th cell
     */
    int getNumberOfPointsOfCell( const int& ima ) const {
      this->_checkCellOOR(ima);
      return this->_ptCell(ima);
    }

    /**
     * @brief Get number of sub-points of the i-th cell
     */
    int getNumberOfSubPointsOfCell( const int& ima ) const {
      this->_checkCellOOR(ima);
      return this->_sptCell(ima);
    }

    /**
     * @brief Get number of components of the i-th cell
     */
    int getNumberOfComponentsForSubpointsOfCell( const int& ima ) const {
      this->_checkCellOOR(ima);
      return this->_cmpsSptCell(ima);
    }

    /**
     * @brief Get number of components
     */
    int getNumberOfComponents() const { return _nbComp; }

    /**
     * @brief Get number of nodes
     */
    int getNumberOfCells() const { return _nbCells; }

    /**
     * @brief Get number of points
     */
    int getMaxNumberOfPoints() const { return _nbPt; }

    /**
     * @brief Get number of sub-points
     */
    int getMaxNumberOfSubPoints() const { return _nbSpt; }

    /**
     * @brief Get the name of the i-th component
     */
    std::string getNameOfComponent( const int& icmp ) const {

      if ( icmp < 0 || icmp >= this->getNumberOfComponents()) {
        throw std::runtime_error( "Component '"+std::to_string(icmp)
                                  +"' is out of range");
      };
      return this->_nameCmp(icmp);
    };

    /**
     * @Brief Get the names of all the components
     */
    VectorString getNameOfComponents() const {

      int size = this->getNumberOfComponents();
      VectorString names;
      names.reserve(size);
      for ( int icmp = 0 ; icmp < size; icmp++ ) {
        names.push_back(this->_nameCmp(icmp));
      }
      return names;
    }

    /**
     * @brief Get physical quantity
     */
    std::string getPhysicalQuantity() const {
      return trim(( *_descriptor )[1].toString());
    }

    /**
     * @brief Get field location
     */
    std::string getFieldLocation() const {
      return trim(( *_descriptor )[2].toString());
    }

    /**
     * @brief Get cells holding components
     */
    VectorLong getCellsWithComponents( ) const {
      VectorLong values;
      for (int ima=0 ; ima<this->getNumberOfCells() ;  ima++){
        if ( this->_cmpsSptCell(ima) > 0 )
          values.push_back(ima);
      }
      return values;
    }

    /**
     * @brief Get values on cells holding components, with mask
     */
    PyObject *getValues( ) {

      int sz = 0;
      int ncmp_max = this->getNumberOfComponents();
      VectorLong cells = this->getCellsWithComponents();

      for ( const auto &ima : cells ) {
        sz = sz + this->_ptCell(ima) * this->_sptCell(ima);
      }
      AS_ASSERT( sz > 0 );

      npy_intp dims[2] = {sz, ncmp_max};
      PyObject *data = PyArray_ZEROS(2, dims, npy_type< ValueType >::value, 0);
      PyObject *mask = PyArray_ZEROS(2, dims, NPY_BOOL, 0);

      ValueType *dataptr = (double *)PyArray_DATA((PyArrayObject *) data);
      bool *maskptr = (bool *)PyArray_DATA((PyArrayObject *) mask);

      int j=0;
      for ( const auto &ima : cells ) {
        for ( int ipt = 0 ; ipt < this->_ptCell(ima); ipt++ ) {
          for ( int ispt = 0 ; ispt < this->_sptCell(ima); ispt++ ) {
            for ( int icmp = 0 ; icmp < this->_cmpsSptCell(ima); icmp++) {
              int posjv = this->_positionInArray(icmp, ima, ipt, ispt);
              if ( ( *_allocated )[posjv] ) {
                dataptr[j+icmp] = ( *_values )[posjv] ;
                maskptr[j+icmp] = true ;
              }
            }
            j = j + ncmp_max;
          }
        }
      }

      PyArray_ENABLEFLAGS((PyArrayObject*) data, NPY_ARRAY_OWNDATA);
      PyArray_ENABLEFLAGS((PyArrayObject*) mask, NPY_ARRAY_OWNDATA);

      PyObject *resu_tuple = PyTuple_New( 2 );
      PyTuple_SetItem(resu_tuple, 0, data);
      PyTuple_SetItem(resu_tuple, 1, mask);

      return resu_tuple;
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
            _nbCells = ( *_size )[0];
            _nbComp = ( *_size )[1];
            _nbPt = ( *_size )[2];
            _nbSpt = ( *_size )[3];
            if ( _values->size() != this->_nbValArray() )
              throw std::runtime_error( "Programming error" );
        }
        return retour;
    };
};

/** @typedef SimpleFieldOnCellsReal Class d'une carte de double */
typedef SimpleFieldOnCells< ASTERDOUBLE > SimpleFieldOnCellsReal;

/**
 * @typedef SimpleFieldOnCellsPtrReal
 * @brief Definition d'un champ aux éléments de double
 */
typedef boost::shared_ptr< SimpleFieldOnCellsReal > SimpleFieldOnCellsRealPtr;

/** @typedef SimpleFieldOnCellsLong Class d'une carte de long */
typedef SimpleFieldOnCells< ASTERINTEGER > SimpleFieldOnCellsLong;

/**
 * @typedef SimpleFieldOnCellsPtrLong
 * @brief Definition d'un champ aux éléments de long
 */
typedef boost::shared_ptr< SimpleFieldOnCellsLong > SimpleFieldOnCellsLongPtr;

#endif /* SIMPLEFIELDONCELLS_H_ */
