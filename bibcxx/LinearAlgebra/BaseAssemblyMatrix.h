#ifndef BASEASSEMBLYMATRIX_H_
#define BASEASSEMBLYMATRIX_H_

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
#include "aster_fort_petsc.h"
#include "astercxx.h"

#include "DataStructures/DataStructure.h"
#include "LinearAlgebra/ElementaryMatrix.h"
#include "Loads/ListOfLoads.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Studies/PhysicalProblem.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Tools.h"

#ifdef ASTER_HAVE_PETSC4PY
#include <petscmat.h>
#endif

#include "Modeling/Model.h"
#include "Numbering/DOFNumbering.h"
#include "Numbering/ParallelDOFNumbering.h"

/**
 * @class BaseAssemblyMatrix
 * @brief Classe template definissant la base d'une sd_matr_asse.
 */

class BaseAssemblyMatrix : public DataStructure {
  protected:
    /** @brief Objet Jeveux '.REFA' */
    JeveuxVectorChar24 _description;
    /** @brief Objet '.CONL' */
    JeveuxVectorReal _scaleFactorLagrangian;
    /** @brief Objet Jeveux '.LIME' */
    JeveuxVectorChar24 _listOfElementaryMatrix;
    /** @brief Objet Jeveux '.PERM' */
    JeveuxVectorLong _perm;

    /** @brief Objet Jeveux '.CCID' */
    JeveuxVectorLong _ccid;
    /** @brief Objet Jeveux '.CCLL' */
    JeveuxVectorLong _ccll;
    /** @brief Objet Jeveux '.CCII' */
    JeveuxVectorLong _ccii;

    /** @brief Objet nume_ddl */
    BaseDOFNumberingPtr _dofNum;
    /** @brief La matrice est elle vide ? */
    bool _isEmpty;
    /** @brief La matrice est elle vide ? */
    bool _isFactorized;
    /** @brief Liste de charges cinematiques */
    ListOfLoadsPtr _listOfLoads;
    /** @brief Solver name (MUMPS or PETSc) */
    std::string _solverName;

  public:
    /**
     * @typedef BaseAssemblyMatrixPtr
     * @brief Pointeur intelligent vers un BaseAssemblyMatrix
     */
    typedef boost::shared_ptr< BaseAssemblyMatrix > BaseAssemblyMatrixPtr;

    /**
     * @brief Constructeur
     */
    BaseAssemblyMatrix() = delete;

    /**
     * @brief Constructeur
     */
    BaseAssemblyMatrix( const std::string &name, const std::string &type );

    /**
     * @brief Constructeur
     */
    BaseAssemblyMatrix( const std::string &type )
        : BaseAssemblyMatrix( ResultNaming::getNewResultName(), type ){};

    /**
     * @brief Constructeur
     */
    BaseAssemblyMatrix( const PhysicalProblemPtr phys_prob, const std::string &type );

    /**
     * @brief Function d'ajout d'un chargement
     * @param Args... Liste d'arguments template
     */
    template < typename... Args > void addLoad( const Args &...a ) {
        _listOfLoads->addLoad( a... );
    };

    /**
     * @brief Assemblage de la matrice
     */
    virtual bool build() { throw std::runtime_error( "Not allowed" ); };

    /**
     * @brief Get the internal DOFNumbering
     * @return Internal DOFNumbering
     */
    BaseDOFNumberingPtr getDOFNumbering() const { return _dofNum; };

    /**
     * @brief Get model
     * @return Internal Model
     */
    ModelPtr getModel() {
        if ( _dofNum ) {
            return _dofNum->getModel();
        }
        return nullptr;
    };

    /**
     * @brief Get mesh
     * @return Internal mesh
     */
    BaseMeshPtr getMesh() {
        if ( _dofNum ) {
            return _dofNum->getMesh();
        }
        return nullptr;
    };

    /**
     * @brief Set new values
     */
    virtual void setValues( const VectorLong idx, const VectorLong jdx, const VectorReal values ) {
        // Template class raises error. It must be specialized in each instanciated class.
        throw std::runtime_error( "Not implemented" );
    };

    /**
     * @brief Transpose
     */
    void transpose() { CALLO_MATR_ASSE_TRANSPOSE( getName() ); };

    /**
     * @brief Transpose and conjugate
     */
    void transposeConjugate() { CALLO_MATR_ASSE_TRANSPOSE_CONJUGATE( getName() ); };

    /**
     * @brief Print the matrix in code_aster format
     */
    void print() const {
        const ASTERINTEGER unit( 6 );
        std::string format( " " );
        CALLO_MATR_ASSE_PRINT( getName(), &unit, format );
    };

    /**
     * @brief Print the matrix in matlab format
     */
    void print( const std::string format ) const {
        const ASTERINTEGER unit( 6 );
        CALLO_MATR_ASSE_PRINT( getName(), &unit, format );
    };

    /**
     * @brief Print the matrix in matlab format in given logical unit
     */
    void print( const ASTERINTEGER unit, const std::string format ) const {
        CALLO_MATR_ASSE_PRINT( getName(), &unit, format );
    };

    /**
     * @brief Get MaterialField
     * @return MaterialField of the first ElementaryMatrix (all others must be the same)
     */
    virtual MaterialFieldPtr getMaterialField() const {
        throw std::runtime_error( "Not allowed" );
    };

    /**
     * @brief Get the number of defined ElementaryMatrix
     * @return size of vector containing ElementaryMatrix
     */
    virtual ASTERINTEGER getNumberOfElementaryMatrix() const {
        throw std::runtime_error( "Not allowed" );
    };

#ifdef ASTER_HAVE_PETSC4PY
    /**
     * @brief Conversion to petsc4py
     * @return converted matrix
     */
    Mat toPetsc4py() {
        Mat myMat;
        PetscErrorCode ierr;

        if ( _isEmpty )
            throw std::runtime_error( "Assembly matrix is empty" );
        if ( getType() != "MATR_ASSE_DEPL_R" )
            throw std::runtime_error( "Not yet implemented" );

        CALLO_MATASS2PETSC( getName(), &myMat, &ierr );

        return myMat;
    };
#endif

    /**
     * @brief Methode permettant de savoir si la matrice est vide
     * @return true si vide
     */
    bool isEmpty() const { return _isEmpty; };

    /**
     * @brief Methode permettant de savoir si la matrice est factorisée
     * @return true si factorisée
     */
    bool isFactorized() const { return _isFactorized; };

    void isFactorized( const bool &facto ) { _isFactorized = facto; };

    /**
     * @brief Methode permettant de definir la numerotation
     * @param currentNum objet DOFNumbering
     */
    void setDOFNumbering( const BaseDOFNumberingPtr currentNum ) { _dofNum = currentNum; };

    /**
     * @brief Function to set the solver name (MUMS or PETSc)
     * @param sName name of solver ("MUMPS" or "PETSC")
     * @todo delete this function and the attribute _solverName
     */
    void setSolverName( const std::string &sName ) { _solverName = sName; };

    /**
     * @brief Delete the factorized matrix used by MUMPS or PETSc if it exist
     * @param sName name of solver ("MUMPS" or "PETSC")
     * @todo delete this function and the attribute _solverName
     */
    bool deleteFactorizedMatrix( void ) {
        if ( _description->exists()  && get_sh_jeveux_status() == 1 ) {
            CALLO_DELETE_MATRIX( getName(), _solverName );
        }

        _isFactorized = false;

        return true;
    };

    /**
     * @brief Return True if CCID object exists for DirichletElimination
     */
    bool hasDirichletEliminationDOFs() const { return _ccid->exists(); }

    /**
     * @brief Return CCID object if exists for DirichletElimination
     */
    JeveuxVectorLong getDirichletBCDOFs() const {
        if ( hasDirichletEliminationDOFs() )
            return _ccid;

        raiseAsterError( "JeveuxError: CCID not existing" );

        return JeveuxVectorLong();
    }

    ASTERDOUBLE
    getLagrangeScaling() const {
        if ( _scaleFactorLagrangian->exists() ) {
            ASTERDOUBLE scaling = 0.0;
            CALLO_CONLAG( getName(), &scaling );
            return scaling;
        }

        return 1.0;
    }

    ListOfLoadsPtr getListOfLoads() const { return _listOfLoads; }

    /**
     * @brief Methode permettant de definir la liste de chargement
     * @param lLoads objet de type ListOfLoadsPtr
     */
    void setListOfLoads( const ListOfLoadsPtr load ) {
        if ( !load )
            raiseAsterError( "Empty load" );

        _listOfLoads = load;
    }

    virtual BaseAssemblyMatrixPtr getEmptyMatrix( const std::string &name ) const {
        throw std::runtime_error( "Not allowed" );
    }
};

typedef boost::shared_ptr< BaseAssemblyMatrix > BaseAssemblyMatrixPtr;

#endif /* BASEASSEMBLYMATRIX_H_ */
