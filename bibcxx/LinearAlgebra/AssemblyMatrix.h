#ifndef ASSEMBLYMATRIX_H_
#define ASSEMBLYMATRIX_H_

/**
 * @file AssemblyMatrix.h
 * @brief Fichier entete de la classe AssemblyMatrix
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

#include "aster_fort_calcul.h"
#include "aster_fort_ds.h"

#include "LinearAlgebra/BaseAssemblyMatrix.h"
#include "LinearAlgebra/ElementaryMatrix.h"
#include "Loads/PhysicalQuantity.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Numbering/DOFNumbering.h"
#include "Numbering/ParallelDOFNumbering.h"
#include "Studies/PhysicalProblem.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Tools.h"

// Forward Declaration
class LinearSolver;
using LinearSolverPtr = std::shared_ptr< LinearSolver >;

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
    typedef std::shared_ptr< ElementaryMatrix< ValueType, PhysicalQuantity > > ElementaryMatrixPtr;

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
    std::map< std::string, ASTERDOUBLE > _coeffReal;

    /** @brief Objet Jeveux '.SOLVEUR' */
    LinearSolverPtr _solver;

    bool isSimilarTo( const AssemblyMatrix &mat ) {
        if ( _dofNum != mat.getDOFNumbering() ) {
            return false;
        }

        return true;
    }

    friend class LinearSolver;

  public:
    /**
     * @typedef AssemblyMatrixPtr
     * @brief Pointeur intelligent vers un AssemblyMatrix
     */
    typedef std::shared_ptr< AssemblyMatrix< ValueType, PhysicalQuantity > > AssemblyMatrixPtr;

    /**
     * @brief Constructeur
     */
    AssemblyMatrix() : AssemblyMatrix( DataStructureNaming::getNewName() ) {};

    /**
     * @brief Constructeur
     */
    AssemblyMatrix( const std::string &name );

    /**
     * @brief Constructeur
     */
    AssemblyMatrix( const PhysicalProblemPtr phys_prob );

    /**
     * @brief Constructeur
     */
    AssemblyMatrix( const std::string &name, const AssemblyMatrix &toCopy );

    AssemblyMatrix( const AssemblyMatrix &toCopy )
        : AssemblyMatrix( DataStructureNaming::getNewName(), toCopy ) {};

    AssemblyMatrix( AssemblyMatrix &&other );

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
    void addElementaryMatrix( const ElementaryMatrixPtr &currentElemMatrix,
                              ASTERDOUBLE coeff = 1.0 ) {
        if ( currentElemMatrix && !currentElemMatrix->isEmpty() ) {
            if ( currentElemMatrix->hasElementaryTerms() ||
                 currentElemMatrix->existsSuperElement() ) {
                _elemMatrix.push_back( currentElemMatrix );
                _coeffReal[currentElemMatrix->getName()] = coeff;
            }
        }
    };

    /**
     * @brief Assemblage de la matrice
     */
    bool assemble( bool clean = true );

    /**
     * @brief Clear all ElementaryMatrixPtr
     */
    void clearElementaryMatrix() {
        _elemMatrix.clear();
        _coeffReal.clear();
    };

    /**
     * @brief Set new values
     */
    void setValues( const VectorLong &idx, const VectorLong &jdx,
                    const std::vector< ValueType > &values ) {
        // Template class raises error. It must be specialized in each instanciated class.
        throw std::runtime_error( "Not implemented" );
    };

    /**
     * @brief Scale the matrix using right and left multiplication
     * by diagonal matrices stored as vectors
     */
    void scale( const std::vector< ValueType > &lvect, const std::vector< ValueType > &rvect ) {
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
        return std::make_shared< AssemblyMatrix< ValueType, PhysicalQuantity > >( name );
    }

    /**
     * @brief Unary Minus overloading
     */
    AssemblyMatrix< ValueType, PhysicalQuantity > operator-() const {
        AssemblyMatrix< ValueType, PhysicalQuantity > tmp( *this );
        tmp *= ASTERDOUBLE( -1 );
        return tmp;
    };

    /**
     * @brief TimesEqual overloading
     * @return Updated field
     */
    AssemblyMatrix< ValueType, PhysicalQuantity > &operator*=( const ASTERDOUBLE &scal ) {

        ( *_matrixValues ) *= scal;

        return *this;
    };

    /**
     * @brief PlusEqual overloading
     * @return Updated field
     */
    AssemblyMatrix< ValueType, PhysicalQuantity > &
    operator+=( const AssemblyMatrix< ValueType, PhysicalQuantity > &mat ) {
        AS_ASSERT( isSimilarTo( mat ) );
        const ASTERDOUBLE c1 = 1.0, c2 = 1.0;
        CALL_ADDMATRASSE( getName(), mat.getName(), &c1, &c2, getName() );

        return *this;
    };

    /**
     * @brief MinusEqual overloading
     * @return Updated field
     */
    AssemblyMatrix< ValueType, PhysicalQuantity > &
    operator-=( const AssemblyMatrix< ValueType, PhysicalQuantity > &mat ) {
        AS_ASSERT( isSimilarTo( mat ) );
        const ASTERDOUBLE c1 = 1.0, c2 = -1.0;
        CALL_ADDMATRASSE( getName(), mat.getName(), &c1, &c2, getName() );

        return *this;
    };

    /**
     * @brief Adding two matrix
     * @return New matrix
     */
    friend AssemblyMatrix< ValueType, PhysicalQuantity >
    operator+( const AssemblyMatrix< ValueType, PhysicalQuantity > &mat1,
               const AssemblyMatrix< ValueType, PhysicalQuantity > &mat2 ) {

        AssemblyMatrix< ValueType, PhysicalQuantity > mat;
        mat.setDOFNumbering( mat1.getDOFNumbering() );
        mat.setListOfLoads( mat1.getListOfLoads() );
        AS_ASSERT( mat.isSimilarTo( mat2 ) );
        const ASTERDOUBLE c1 = 1.0, c2 = 1.0;
        CALL_ADDMATRASSE( mat1.getName(), mat2.getName(), &c1, &c2, mat.getName() );
        mat._isEmpty = false;
        return mat;
    };

    /**
     * @brief Adding two matrix
     * @return New matrix
     */
    friend AssemblyMatrix< ValueType, PhysicalQuantity >
    operator-( const AssemblyMatrix< ValueType, PhysicalQuantity > &mat1,
               const AssemblyMatrix< ValueType, PhysicalQuantity > &mat2 ) {
        AssemblyMatrix< ValueType, PhysicalQuantity > mat;
        mat.setDOFNumbering( mat1.getDOFNumbering() );
        mat.setListOfLoads( mat1.getListOfLoads() );
        AS_ASSERT( mat.isSimilarTo( mat2 ) );
        const ASTERDOUBLE c1 = 1.0, c2 = -1.0;
        CALL_ADDMATRASSE( mat1.getName(), mat2.getName(), &c1, &c2, mat.getName() );
        mat._isEmpty = false;
        return mat;
    };

    /**
     * @brief Multiply by a scalar on left overloading
     * @return New matrix
     */
    friend AssemblyMatrix< ValueType, PhysicalQuantity >
    operator*( const ASTERDOUBLE &scal, AssemblyMatrix< ValueType, PhysicalQuantity > mat ) {

        mat *= scal;

        return mat;
    };

    /**
     * @brief Multiply by a matrix by a vector
     * @return New field
     */
    friend FieldOnNodes< ValueType >
    operator*( const AssemblyMatrix< ValueType, PhysicalQuantity > &mat,
               const FieldOnNodes< ValueType > &lhs ) {

        auto matDesc = mat.getDOFNumbering()->getEquationNumbering();
        auto vecDesc = lhs.getDescription();
        if ( matDesc != vecDesc ) {
            if ( *matDesc != *vecDesc ) {
                AS_ABORT( "Incompatible numbering: " + matDesc->getName() + " vs " +
                          vecDesc->getName() );
            }
        }
        FieldOnNodes< ValueType > result( mat.getDOFNumbering() );

        CALL_MVMULT( mat.getName(), lhs.getName(), result.getName() );

        return result;
    };

    /**
     * @brief Apply the DirichletBC into the Rhs (aka kinematic aka no Lagrange multipliers)
     * @param DirichletBC [FieldOnNodes] the values on the DirichletBC
     * @param Rhs [FieldOnNodes] The residual to be modified
     */
    void applyDirichletBC( const FieldOnNodes< ValueType > &DirichletBC,
                           FieldOnNodes< ValueType > &Rhs ) const {
        // Template class raises error. It must be specialized in each instanciated class.
        throw std::runtime_error( "Not implemented" );
    };

    AssemblyMatrix copy() { return *this; };

    /**
     * @brief Methode permettant de definir si un solveur est attribué à la matrice
     */
    void defineSolver();

    std::vector< ValueType > getUpperValues() const {
        _matrixValues->build();
        _matrixValues->updateValuePointer();
        return ( *_matrixValues )[1]->toVector();
    };

    std::vector< ValueType > getLowerValues() const {
        _matrixValues->build();
        _matrixValues->updateValuePointer();
        return ( *_matrixValues )[2]->toVector();
    };
};

/** @typedef Definition d'une matrice assemblee de double */
template <>
void AssemblyMatrix< ASTERDOUBLE, Displacement >::setValues( const VectorLong &idx,
                                                             const VectorLong &jdx,
                                                             const VectorReal &values );
template <>
void AssemblyMatrix< ASTERDOUBLE, Displacement >::scale( const VectorReal &lvect,
                                                         const VectorReal &rvect );
template <>
void AssemblyMatrix< ASTERDOUBLE, Displacement >::applyDirichletBC(
    const FieldOnNodesReal &DirichletBC, FieldOnNodesReal &Rhs ) const;

typedef AssemblyMatrix< ASTERDOUBLE, Displacement > AssemblyMatrixDisplacementReal;

/** @typedef Definition d'une matrice assemblee de double avec relations linéaires éliminées*/
class AssemblyMatrixEliminatedReal : public AssemblyMatrixDisplacementReal {
  public:
    using AssemblyMatrixDisplacementReal::AssemblyMatrixDisplacementReal;
};

/** @typedef Definition d'une matrice assemblee de complexe */
template class AssemblyMatrix< ASTERCOMPLEX, Displacement >;
typedef AssemblyMatrix< ASTERCOMPLEX, Displacement > AssemblyMatrixDisplacementComplex;

/** @typedef Definition d'une matrice assemblee de double temperature */
template <>
void AssemblyMatrix< ASTERDOUBLE, Temperature >::setValues( const VectorLong &idx,
                                                            const VectorLong &jdx,
                                                            const VectorReal &values );
typedef AssemblyMatrix< ASTERDOUBLE, Temperature > AssemblyMatrixTemperatureReal;

template <>
void AssemblyMatrix< ASTERDOUBLE, Temperature >::applyDirichletBC(
    const FieldOnNodesReal &DirichletBC, FieldOnNodesReal &Rhs ) const;

/** @typedef Definition d'une matrice assemblee de double pression */
template <>
void AssemblyMatrix< ASTERDOUBLE, Pressure >::setValues( const VectorLong &idx,
                                                         const VectorLong &jdx,
                                                         const VectorReal &values );
typedef AssemblyMatrix< ASTERDOUBLE, Pressure > AssemblyMatrixPressureReal;

template <>
void AssemblyMatrix< ASTERDOUBLE, Pressure >::applyDirichletBC( const FieldOnNodesReal &DirichletBC,
                                                                FieldOnNodesReal &Rhs ) const;

/** @typedef Definition d'une matrice assemblee de ASTERCOMPLEX temperature */
template class AssemblyMatrix< ASTERCOMPLEX, Temperature >;
typedef AssemblyMatrix< ASTERCOMPLEX, Temperature > AssemblyMatrixTemperatureComplex;

/** @typedef Definition d'une matrice assemblee de ASTERCOMPLEX pression */
template class AssemblyMatrix< ASTERCOMPLEX, Pressure >;
typedef AssemblyMatrix< ASTERCOMPLEX, Pressure > AssemblyMatrixPressureComplex;

typedef std::shared_ptr< AssemblyMatrixDisplacementReal > AssemblyMatrixDisplacementRealPtr;
typedef std::shared_ptr< AssemblyMatrixDisplacementComplex > AssemblyMatrixDisplacementComplexPtr;
typedef std::shared_ptr< AssemblyMatrixTemperatureReal > AssemblyMatrixTemperatureRealPtr;
typedef std::shared_ptr< AssemblyMatrixTemperatureComplex > AssemblyMatrixTemperatureComplexPtr;
typedef std::shared_ptr< AssemblyMatrixPressureReal > AssemblyMatrixPressureRealPtr;
typedef std::shared_ptr< AssemblyMatrixPressureComplex > AssemblyMatrixPressureComplexPtr;
typedef std::shared_ptr< AssemblyMatrixEliminatedReal > AssemblyMatrixEliminatedRealPtr;

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
      _ccva( JeveuxVector< ValueType >( getName() + ".CCVA" ) ) {};

template < class ValueType, PhysicalQuantityEnum PhysicalQuantity >
AssemblyMatrix< ValueType, PhysicalQuantity >::AssemblyMatrix( const PhysicalProblemPtr phys_prob )
    : AssemblyMatrix() {
    _dofNum = phys_prob->getDOFNumbering();
    _listOfLoads = phys_prob->getListOfLoads();
};

template < class ValueType, PhysicalQuantityEnum PhysicalQuantity >
AssemblyMatrix< ValueType, PhysicalQuantity >::AssemblyMatrix( const std::string &name,
                                                               const AssemblyMatrix &toCopy )
    : BaseAssemblyMatrix( name,
                          "MATR_ASSE_" + std::string( PhysicalQuantityNames[PhysicalQuantity] ) +
                              ( typeid( ValueType ) == typeid( ASTERDOUBLE ) ? "_R" : "_C" ),
                          toCopy ),
      _matrixValues( JeveuxCollection< ValueType >( getName() + ".VALM" ) ),
      _valf( JeveuxVector< ValueType >( getName() + ".VALF" ) ),
      _walf( JeveuxVector< ValueType >( getName() + ".WALF" ) ),
      _ualf( JeveuxVector< ValueType >( getName() + ".UALF" ) ),
      _digs( JeveuxVector< ValueType >( getName() + ".DIGS" ) ),
      _ccva( JeveuxVector< ValueType >( getName() + ".CCVA" ) ) {

    ( *_matrixValues ) = ( *toCopy._matrixValues );
    ( *_valf ) = ( *toCopy._valf );
    ( *_walf ) = ( *toCopy._walf );
    ( *_ualf ) = ( *toCopy._ualf );
    ( *_digs ) = ( *toCopy._digs );
    ( *_ccva ) = ( *toCopy._ccva );

    _elemMatrix = toCopy._elemMatrix;
    _coeffReal = toCopy._coeffReal;
}

template < class ValueType, PhysicalQuantityEnum PhysicalQuantity >
AssemblyMatrix< ValueType, PhysicalQuantity >::AssemblyMatrix( AssemblyMatrix &&other )
    : BaseAssemblyMatrix( std::move( other ) ) {

    _matrixValues = other._matrixValues;
    _valf = other._valf;
    _walf = other._walf;
    _ualf = other._ualf;
    _digs = other._digs;
    _ccva = other._ccva;

    _elemMatrix = std::move( other._elemMatrix );
    _coeffReal = std::move( other._coeffReal );
}

template < class ValueType, PhysicalQuantityEnum PhysicalQuantity >
bool AssemblyMatrix< ValueType, PhysicalQuantity >::assemble( bool clean ) {
    if ( !_dofNum->exists() )
        throw std::runtime_error( "Numbering is empty" );

    if ( getNumberOfElementaryMatrix() == 0 )
        throw std::runtime_error( "Elementary matrix is empty" );

    ASTERINTEGER typscal = 2;
    if ( typeid( ValueType ) == typeid( ASTERDOUBLE ) )
        typscal = 1;

    JeveuxVectorReal list_coef( "&&AM.COEF_R" );

    ASTERINTEGER nbMatrElem = _elemMatrix.size();
    VectorString names;
    names.reserve( nbMatrElem );
    list_coef->reserve( nbMatrElem );
    for ( const auto elemIt : _elemMatrix ) {
        auto matr = elemIt->getName();
        names.push_back( matr );
        list_coef->push_back( _coeffReal[matr] );
    }

    char *tabNames = vectorStringAsFStrArray( names, 19 );

    std::string base( "G" );
    std::string cumul( "ZERO" );

    if ( _listOfLoads->isEmpty() && _listOfLoads->getNumberOfLoads() != 0 )
        _listOfLoads->build();

    CALL_ASMATR( &nbMatrElem, tabNames, list_coef->getName().c_str(), _dofNum->getName().c_str(),
                 _listOfLoads->getName().c_str(), cumul.c_str(), base.c_str(), &typscal,
                 getName().c_str() );
    _isEmpty = false;

    // free matr_elem string
    FreeStr( tabNames );

    if ( !isMPIFull() && !getMesh()->isParallel() ) {
        std::string type = "MATR_ASSE";
        CALLO_SDMPIC( type, getName() );
    }

    if ( clean ) {
        clearElementaryMatrix();
    }

    return true;
};

#endif /* ASSEMBLYMATRIX_H_ */
