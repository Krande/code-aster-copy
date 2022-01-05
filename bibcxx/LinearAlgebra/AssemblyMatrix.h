#ifndef ASSEMBLYMATRIX_H_
#define ASSEMBLYMATRIX_H_

/**
 * @file AssemblyMatrix.h
 * @brief Fichier entete de la classe AssemblyMatrix
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

#include "aster_fort_calcul.h"
#include "aster_fort_ds.h"
#include "astercxx.h"

#include "LinearAlgebra/ElementaryMatrix.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Studies/PhysicalProblem.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Tools.h"

#include "LinearAlgebra/BaseAssemblyMatrix.h"
#include "Loads/PhysicalQuantity.h"
#include "Modeling/Model.h"
#include "Numbering/DOFNumbering.h"
#include "Numbering/ParallelDOFNumbering.h"

/**
 * @class AssemblyMatrix
 * @brief Classe template definissant une sd_matr_asse.
 *        Cette classe est volontairement succinte car on n'en connait pas encore l'usage
 * @author Nicolas Sellenet
 * @todo revoir le template pour prendre la grandeur en plus
 */
template < class ValueType, PhysicalQuantityEnum PhysicalQuantity >
class AssemblyMatrix : public BaseAssemblyMatrix {
  private:
    typedef boost::shared_ptr< ElementaryMatrix< ValueType, PhysicalQuantity > >
        ElementaryMatrixPtr;

    /** @brief Collection '.VALM' */
    JeveuxCollection< ValueType > _matrixValues;
    /** @brief Objet Jeveux '.VALF' */
    JeveuxVector< ValueType > _valf;
    /** @brief Objet Jeveux '.WALF' */
    JeveuxVector< ValueType > _walf;
    /** @brief Objet Jeveux '.UALF' */
    JeveuxVector< ValueType > _ualf;
    /** @brief Objet Jeveux '.DIGS' */
    JeveuxVector< ValueType > _digs;
    /** @brief Objet Jeveux '.CCVA' */
    JeveuxVector< ValueType > _ccva;

    /** @brief ElementaryMatrix sur lesquelles sera construit la matrice */
    std::vector< ElementaryMatrixPtr > _elemMatrix;

  public:
    /**
     * @typedef AssemblyMatrixPtr
     * @brief Pointeur intelligent vers un AssemblyMatrix
     */
    typedef boost::shared_ptr< AssemblyMatrix< ValueType, PhysicalQuantity > > AssemblyMatrixPtr;

    /**
     * @brief Constructeur
     */
    AssemblyMatrix() : AssemblyMatrix( ResultNaming::getNewResultName() ){};

    /**
     * @brief Constructeur
     */
    AssemblyMatrix( const std::string &name );

    /**
     * @brief Constructeur
     */
    AssemblyMatrix( const PhysicalProblemPtr phys_prob );

    /**
     * @brief Destructeur
     */
    ~AssemblyMatrix() {
        // #ifdef ASTER_DEBUG_CXX
        //         std::cout << "DEBUG: BaseAssemblyMatrix.destr: " << this->getName() << std::endl;
        // #endif
        // two temporary objects to delete
        CALLO_JEDETR( getName() + ".&INT" );
        CALLO_JEDETR( getName() + ".&IN2" );
        this->deleteFactorizedMatrix();
    };

    /**
     * @brief Methode permettant de definir les matrices elementaires
     * @param currentElemMatrix objet ElementaryMatrix
     */
    void appendElementaryMatrix( const ElementaryMatrixPtr &currentElemMatrix ) {
        _elemMatrix.push_back( currentElemMatrix );
    };

    /**
     * @brief Assemblage de la matrice
     */
    bool build();

    /**
     * @brief Clear all ElementaryMatrixPtr
     */
    void clearElementaryMatrix() { _elemMatrix.clear(); };

    /**
     * @brief Set new values
     */
    void setValues( const VectorLong idx, const VectorLong jdx,
                    const std::vector< ValueType > values ) {
        // Template class raises error. It must be specialized in each instanciated class.
        throw std::runtime_error( "Not implemented" );
    };

    /**
     * @brief Get MaterialField
     * @return MaterialField of the first ElementaryMatrix (all others must be the same)
     */
    MaterialFieldPtr getMaterialField() const {
        if ( _elemMatrix.size() != 0 )
            return _elemMatrix[0]->getMaterialField();
        throw std::runtime_error( "No ElementaryMatrix in AssemblyMatrix" );
    };

    /**
     * @brief Get the number of defined ElementaryMatrix
     * @return size of vector containing ElementaryMatrix
     */
    ASTERINTEGER getNumberOfElementaryMatrix() const { return _elemMatrix.size(); };

    BaseAssemblyMatrixPtr getEmptyMatrix( const std::string &name ) const {
        return boost::make_shared< AssemblyMatrix< ValueType, PhysicalQuantity > >( name );
    }
};

/** @typedef Definition d'une matrice assemblee de double */
template <>
void AssemblyMatrix< ASTERDOUBLE, Displacement >::setValues( const VectorLong idx,
                                                             const VectorLong jdx,
                                                             const VectorReal values );
typedef AssemblyMatrix< ASTERDOUBLE, Displacement > AssemblyMatrixDisplacementReal;

/** @typedef Definition d'une matrice assemblee de complexe */
template class AssemblyMatrix< ASTERCOMPLEX, Displacement >;
typedef AssemblyMatrix< ASTERCOMPLEX, Displacement > AssemblyMatrixDisplacementComplex;

/** @typedef Definition d'une matrice assemblee de double temperature */
template <>
void AssemblyMatrix< ASTERDOUBLE, Temperature >::setValues( const VectorLong idx,
                                                            const VectorLong jdx,
                                                            const VectorReal values );
typedef AssemblyMatrix< ASTERDOUBLE, Temperature > AssemblyMatrixTemperatureReal;

/** @typedef Definition d'une matrice assemblee de double pression */
template <>
void AssemblyMatrix< ASTERDOUBLE, Pressure >::setValues( const VectorLong idx, const VectorLong jdx,
                                                         const VectorReal values );
typedef AssemblyMatrix< ASTERDOUBLE, Pressure > AssemblyMatrixPressureReal;

/** @typedef Definition d'une matrice assemblee de ASTERCOMPLEX temperature */
template class AssemblyMatrix< ASTERCOMPLEX, Temperature >;
typedef AssemblyMatrix< ASTERCOMPLEX, Temperature > AssemblyMatrixTemperatureComplex;

/** @typedef Definition d'une matrice assemblee de ASTERCOMPLEX pression */
template class AssemblyMatrix< ASTERCOMPLEX, Pressure >;
typedef AssemblyMatrix< ASTERCOMPLEX, Pressure > AssemblyMatrixPressureComplex;

typedef boost::shared_ptr< AssemblyMatrixDisplacementReal > AssemblyMatrixDisplacementRealPtr;
typedef boost::shared_ptr< AssemblyMatrixDisplacementComplex > AssemblyMatrixDisplacementComplexPtr;
typedef boost::shared_ptr< AssemblyMatrixTemperatureReal > AssemblyMatrixTemperatureRealPtr;
typedef boost::shared_ptr< AssemblyMatrixTemperatureComplex > AssemblyMatrixTemperatureComplexPtr;
typedef boost::shared_ptr< AssemblyMatrixPressureReal > AssemblyMatrixPressureRealPtr;
typedef boost::shared_ptr< AssemblyMatrixPressureComplex > AssemblyMatrixPressureComplexPtr;

template < class ValueType, PhysicalQuantityEnum PhysicalQuantity >
AssemblyMatrix< ValueType, PhysicalQuantity >::AssemblyMatrix( const std::string &name )
    : BaseAssemblyMatrix( name,
                          "MATR_ASSE_" + std::string( PhysicalQuantityNames[PhysicalQuantity] ) +
                              ( typeid( ValueType ) == typeid( ASTERDOUBLE ) ? "_R" : "_C" ) ),
      _matrixValues( JeveuxCollection< ValueType >( getName() + ".VALM" ) ),
      _valf( JeveuxVector< ValueType >( getName() + ".VALF" ) ),
      _walf( JeveuxVector< ValueType >( getName() + ".WALF" ) ),
      _ualf( JeveuxVector< ValueType >( getName() + ".UALF" ) ),
      _digs( JeveuxVector< ValueType >( getName() + ".DIGS" ) ),
      _ccva( JeveuxVector< ValueType >( getName() + ".CCVA" ) ){};

template < class ValueType, PhysicalQuantityEnum PhysicalQuantity >
AssemblyMatrix< ValueType, PhysicalQuantity >::AssemblyMatrix( const PhysicalProblemPtr phys_prob )
    : AssemblyMatrix() {
    _dofNum = phys_prob->getDOFNumbering();
    _listOfLoads = phys_prob->getListOfLoads();
};

template < class ValueType, PhysicalQuantityEnum PhysicalQuantity >
bool AssemblyMatrix< ValueType, PhysicalQuantity >::build() {
    if ( _dofNum->isEmpty() )
        throw std::runtime_error( "Numbering is empty" );

    if ( getNumberOfElementaryMatrix() == 0 )
        throw std::runtime_error( "Elementary matrix is empty" );

    ASTERINTEGER typscal = 2;
    if ( typeid( ValueType ) == typeid( ASTERDOUBLE ) )
        typscal = 1;

    ASTERINTEGER nbMatrElem = _elemMatrix.size();
    VectorString names;
    names.reserve( nbMatrElem );
    for ( const auto elemIt : _elemMatrix )
        names.push_back( elemIt->getName() );

    char *tabNames = vectorStringAsFStrArray( names, 8 );

    std::string base( "G" );
    std::string blanc( " " );
    std::string cumul( "ZERO" );

    _listOfLoads->isEmpty();
    if ( _listOfLoads->isEmpty() && _listOfLoads->getNumberOfLoads() != 0 )
        _listOfLoads->build();

    CALL_ASMATR( &nbMatrElem, tabNames, blanc.c_str(), _dofNum->getName().c_str(),
                 _listOfLoads->getName().c_str(), cumul.c_str(), base.c_str(), &typscal,
                 getName().c_str() );
    _isEmpty = false;

    // free matr_elem string
    FreeStr( tabNames );

    return true;
};

#endif /* ASSEMBLYMATRIX_H_ */
