#ifndef LINEARSOLVER_H_
#define LINEARSOLVER_H_

/**
 * @file LinearSolver.h
 * @brief Fichier entete de la classe LinearSolver
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

#include <list>
#include <set>
#include <stdexcept>
#include <string>

#include "DataFields/FieldOnNodes.h"
#include "DataStructures/DataStructure.h"
#include "LinearAlgebra/AssemblyMatrix.h"
#include "MemoryManager/JeveuxVector.h"
#include "Supervis/ResultNaming.h"
#include "astercxx.h"

/**
 * @class LinearSolver
 * @brief Cette classe permet de definir un solveur lineaire
 */
class LinearSolver : public DataStructure {
  protected:
    bool _isEmpty;

    JeveuxVectorChar24 _charValues;
    JeveuxVectorReal _doubleValues;
    JeveuxVectorLong _integerValues;
    JeveuxVectorChar80 _petscOptions;

    AssemblyMatrixDisplacementRealPtr _matrixPrec;
    std::string _commandName;
    bool _xfem;
    PyObject *_keywords = NULL;

  public:
    /**
     * @typedef LinearSolverPtr
     * @brief Pointeur intelligent vers un LinearSolver
     */
    typedef boost::shared_ptr< LinearSolver > LinearSolverPtr;

    /**
     * @brief Constructeur
     */
    LinearSolver() : LinearSolver( ResultNaming::getNewResultName() ){};

    /**
     * @brief Constructeur
     * @param name Name of the DataStructure
     */
    LinearSolver( const std::string name );

    /**
     * @brief Destructor
     */
    ~LinearSolver() { Py_XDECREF( _keywords ); };

    /**
     * @brief Return the solver name.
     * @return string
     */
    // can not be pure virtual because of boost wrapping
    virtual const std::string getSolverName() const { return ""; };

    /**
     * @brief Tell if the solver support HPC distributed parallelism.
     * @return bool
     */
    virtual const bool supportParallelMesh() const { return false; };

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
    void setKeywords( PyObject *user_keywords );

    /**
     * @brief Returns a dict containing the SOLVEUR keyword.
     * @return PyDict (new reference)
     */
    PyObject *getKeywords() const;

    /**
     * @brief Methode permettant de savoir si la matrice est vide
     * @return true si vide
     */
    bool isEmpty() { return _isEmpty; };

    /**
     * @brief Factorisation d'une matrice
     * @param currentMatrix Matrice assemblee
     */
    bool factorize( AssemblyMatrixDisplacementRealPtr currentMatrix );

    /**
     * @brief Inversion du systeme lineaire
     * @param currentMatrix Matrice assemblee
     * @param dirichletBCField Charge cinématique
     * @param currentRHS Second membre
     * @param result champ aux noeuds résultat (optionnel)
     * @return champ aux noeuds resultat
     */
    FieldOnNodesRealPtr
    solve( const AssemblyMatrixDisplacementRealPtr &currentMatrix,
           const FieldOnNodesRealPtr &currentRHS,
           FieldOnNodesRealPtr result = FieldOnNodesRealPtr( new FieldOnNodesReal() ) ) const;

    /**
     * @brief Inversion du systeme lineaire
     * @param currentMatrix Matrice assemblee
     * @param dirichletBCField Charge cinématique
     * @param currentRHS Second membre
     * @param result champ aux noeuds résultat (optionnel)
     * @return champ aux noeuds resultat
     */
    FieldOnNodesRealPtr solveWithDirichletBC(
        const AssemblyMatrixDisplacementRealPtr &currentMatrix,
        const FieldOnNodesRealPtr &dirichletBCField, const FieldOnNodesRealPtr &currentRHS,
        FieldOnNodesRealPtr result = FieldOnNodesRealPtr( new FieldOnNodesReal() ) ) const;

    friend class LinearStaticAnalysis;
};

/**
 * @typedef LinearSolverPtr
 * @brief Pointeur intelligent vers un LinearSolver
 */
typedef boost::shared_ptr< LinearSolver > LinearSolverPtr;

class LdltSolver : public LinearSolver {
  public:
    LdltSolver( const std::string name ) : LinearSolver( name ){};
    LdltSolver() : LinearSolver(){};
    const std::string getSolverName() const { return "LDLT"; };
};

class MultFrontSolver : public LinearSolver {
  public:
    MultFrontSolver( const std::string name ) : LinearSolver( name ){};
    MultFrontSolver() : LinearSolver(){};

    const std::string getSolverName() const { return "MULT_FRONT"; };
};

class MumpsSolver : public LinearSolver {
  public:
    MumpsSolver( const std::string name ) : LinearSolver( name ){};
    MumpsSolver() : LinearSolver(){};
    const std::string getSolverName() const { return "MUMPS"; };
    const bool supportParallelMesh() const { return true; };
};

class PetscSolver : public LinearSolver {
  public:
    PetscSolver( const std::string name ) : LinearSolver( name ){};
    PetscSolver() : LinearSolver(){};
    const std::string getSolverName() const { return "PETSC"; };
    const bool supportParallelMesh() const { return true; };
};

class GcpcSolver : public LinearSolver {
  public:
    GcpcSolver( const std::string name ) : LinearSolver( name ){};
    GcpcSolver() : LinearSolver(){};
    const std::string getSolverName() const { return "GCPC"; };
};

/** @brief Enveloppe d'un pointeur intelligent vers un LinearSolver< MultFront > */
typedef boost::shared_ptr< MultFrontSolver > MultFrontSolverPtr;

/** @brief Enveloppe d'un pointeur intelligent vers un LinearSolver< Ldlt > */
typedef boost::shared_ptr< LdltSolver > LdltSolverPtr;

/** @brief Enveloppe d'un pointeur intelligent vers un LinearSolver< Mumps > */
typedef boost::shared_ptr< MumpsSolver > MumpsSolverPtr;

/** @brief Enveloppe d'un pointeur intelligent vers un LinearSolver< Petsc > */
typedef boost::shared_ptr< PetscSolver > PetscSolverPtr;

/** @brief Enveloppe d'un pointeur intelligent vers un LinearSolver< Gcpc > */
typedef boost::shared_ptr< GcpcSolver > GcpcSolverPtr;

#endif /* LINEARSOLVER_H_ */
