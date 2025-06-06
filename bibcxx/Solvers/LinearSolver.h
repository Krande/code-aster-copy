#ifndef LINEARSOLVER_H_
#define LINEARSOLVER_H_

/**
 * @file LinearSolver.h
 * @brief Fichier entete de la classe LinearSolver
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

#include "astercxx.h"

#include "DataFields/FieldOnNodes.h"
#include "DataStructures/DataStructure.h"
#include "LinearAlgebra/BaseAssemblyMatrix.h"
#include "MemoryManager/JeveuxVector.h"
#include "Supervis/ResultNaming.h"

/**
 * @class LinearSolver
 * @brief Cette classe permet de definir un solveur lineaire
 */
class LinearSolver : public DataStructure {
  protected:
    bool _isBuilt;

    JeveuxVectorChar24 _charValues;
    JeveuxVectorReal _doubleValues;
    JeveuxVectorLong _integerValues;
    JeveuxVectorChar80 _petscOptions;

    BaseAssemblyMatrixPtr _matrix;
    BaseAssemblyMatrixPtr _matrixPrec;
    std::string _cataPath;
    bool _xfem;
    py::object _keywords;

    void _solve( const std::string &rhsName, const std::string &diriName,
                 const std::string &resultName ) const;

  public:
    /**
     * @typedef LinearSolverPtr
     * @brief Pointeur intelligent vers un LinearSolver
     */
    typedef std::shared_ptr< LinearSolver > LinearSolverPtr;

    /**
     * @brief Constructeur
     */
    LinearSolver() : LinearSolver( ResultNaming::getNewResultName() ) {};

    /**
     * @brief Constructeur
     * @param name Name of the DataStructure
     */
    LinearSolver( const std::string name );

    /**
     * @brief Destructor
     */
    ~LinearSolver() { deleteFactorizedMatrix(); };

    /**
     * @brief Return the solver name.
     * @return string
     */
    // can not be pure virtual because of Python wrapping
    virtual std::string getSolverName() const { return ""; };

    /**
     * @brief Tell if the solver support HPC distributed parallelism.
     * @return bool
     */
    virtual bool supportParallelMesh() const { return false; };

    /**
     * @brief Construction de la sd_solveur
     * @return vrai si tout s'est bien passé
     */
    bool build();

    /**
     * @brief Enable Xfem preconditioning
     */
    void enableXfem() { _xfem = true; };

    /**
     * @brief Store user keywords for SOLVEUR.
     */
    void setKeywords( py::object &user_keywords );

    /**
     * @brief Set command name.
     */
    void setCataPath( const std::string &cataPath ) { _cataPath = cataPath; };

    /**
     * @brief Returns a dict containing the SOLVEUR keyword.
     * @return PyDict (new reference)
     */
    py::dict getKeywords() const;

    /**
     * @brief Methode permettant de savoir si la matrice a été construite
     * @return true si la matrice a été construite
     */
    bool isBuilt() { return _isBuilt; };

    /**
     * @brief Factorisation d'une matrice
     * @param currentMatrix Matrice assemblee
     * @return bool True if ok
     */
    bool factorize( const BaseAssemblyMatrixPtr currentMatrix, bool raiseException = false );

    /**
     * @brief Factorisation d'une matrice
     */
    bool deleteFactorizedMatrix();

    /**
     * @brief Get factorized matrix
     */
    BaseAssemblyMatrixPtr getMatrix() const { return _matrix; };

    /**
     * @brief Get Preconditionning matrix
     */
    BaseAssemblyMatrixPtr getPrecondMatrix() const { return _matrixPrec; };

    /**
     * @brief Inversion du systeme lineaire
     * @param currentRHS Second membre
     * @param dirichletBCField Charge cinématique
     * @return champ aux noeuds resultat
     */
    FieldOnNodesRealPtr solve( const FieldOnNodesRealPtr currentRHS,
                               const FieldOnNodesRealPtr dirichletBCField = nullptr ) const;

    FieldOnNodesComplexPtr solve( const FieldOnNodesComplexPtr currentRHS,
                                  const FieldOnNodesComplexPtr dirichletBCField = nullptr ) const;
};

/**
 * @typedef LinearSolverPtr
 * @brief Pointeur intelligent vers un LinearSolver
 */
typedef std::shared_ptr< LinearSolver > LinearSolverPtr;

class LdltSolver : public LinearSolver {
  public:
    LdltSolver( const std::string name ) : LinearSolver( name ) {};
    LdltSolver() : LinearSolver() {};
    std::string getSolverName() const { return "LDLT"; };
};

class MultFrontSolver : public LinearSolver {
  public:
    MultFrontSolver( const std::string name ) : LinearSolver( name ) {};
    MultFrontSolver() : LinearSolver() {};

    std::string getSolverName() const { return "MULT_FRONT"; };
};

class MumpsSolver : public LinearSolver {
  public:
    MumpsSolver( const std::string name ) : LinearSolver( name ) {};
    MumpsSolver() : LinearSolver() {};
    std::string getSolverName() const { return "MUMPS"; };
    bool supportParallelMesh() const { return true; };
};

class PetscSolver : public LinearSolver {
  public:
    PetscSolver( const std::string name ) : LinearSolver( name ) {};
    PetscSolver() : LinearSolver() {};
    std::string getSolverName() const { return "PETSC"; };
    bool supportParallelMesh() const { return true; };
    /**
     * @brief Return the PETSc options set by the user
     * @return string
     */
    std::string getPetscOptions() const {
        std::string opt;
        _petscOptions->updateValuePointer();
        for ( auto i = 0; i < _petscOptions->size(); i++ ) {
            opt.append( ( *_petscOptions )[i] );
        }
        return strip( opt );
    };
};

class GcpcSolver : public LinearSolver {
  public:
    GcpcSolver( const std::string name ) : LinearSolver( name ) {};
    GcpcSolver() : LinearSolver() {};
    std::string getSolverName() const { return "GCPC"; };
};

/** @brief Enveloppe d'un pointeur intelligent vers un LinearSolver< MultFront > */
typedef std::shared_ptr< MultFrontSolver > MultFrontSolverPtr;

/** @brief Enveloppe d'un pointeur intelligent vers un LinearSolver< Ldlt > */
typedef std::shared_ptr< LdltSolver > LdltSolverPtr;

/** @brief Enveloppe d'un pointeur intelligent vers un LinearSolver< Mumps > */
typedef std::shared_ptr< MumpsSolver > MumpsSolverPtr;

/** @brief Enveloppe d'un pointeur intelligent vers un LinearSolver< Petsc > */
typedef std::shared_ptr< PetscSolver > PetscSolverPtr;

/** @brief Enveloppe d'un pointeur intelligent vers un LinearSolver< Gcpc > */
typedef std::shared_ptr< GcpcSolver > GcpcSolverPtr;

#endif /* LINEARSOLVER_H_ */
