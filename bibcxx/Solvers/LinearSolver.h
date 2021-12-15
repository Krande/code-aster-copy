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
#include "Solvers/AllowedLinearSolver.h"
#include "Supervis/ResultNaming.h"
#include "Utilities/GenericParameter.h"
#include "astercxx.h"

/**
 * @struct LinearStaticAnalysis predefinition
 * @todo to remove
 */
class LinearStaticAnalysis;

/* person_in_charge: nicolas.sellenet at edf.fr */

// Ces wrappers sont la pour autoriser que les set soient const
// Sinon, on aurait pas pu passer directement des const set<> en parametre template
/**
 * @struct WrapMultFront
 * @brief Structure destinee a contenir les renumeroteurs autorises pour MultFront
 */
struct WrapMultFront {
    static const LinearSolverEnum solverType = MultFront;
    static const std::set< Renumbering > setOfAllowedRenumbering;
    static const bool isHPCCompliant = false;
    static const Renumbering defaultRenumbering = Metis;
};

/**
 * @struct WrapLdlt
 * @brief Structure destinee a contenir les renumeroteurs autorises pour Ldlt
 */
struct WrapLdlt {
    static const LinearSolverEnum solverType = Ldlt;
    static const std::set< Renumbering > setOfAllowedRenumbering;
    static const bool isHPCCompliant = false;
    static const Renumbering defaultRenumbering = RCMK;
};

/**
 * @struct WrapMumps
 * @brief Structure destinee a contenir les renumeroteurs autorises pour Mumps
 */
struct WrapMumps {
    static const LinearSolverEnum solverType = Mumps;
    static const std::set< Renumbering > setOfAllowedRenumbering;
    static const bool isHPCCompliant = true;
    static const Renumbering defaultRenumbering = Auto;
};

/**
 * @struct WrapPetsc
 * @brief Structure destinee a contenir les renumeroteurs autorises pour Petsc
 */
struct WrapPetsc {
    static const LinearSolverEnum solverType = Petsc;
    static const std::set< Renumbering > setOfAllowedRenumbering;
    static const bool isHPCCompliant = true;
    static const Renumbering defaultRenumbering = Sans;
};

/**
 * @struct WrapGcpc
 * @brief Structure destinee a contenir les renumeroteurs autorises pour Gcpc
 */
struct WrapGcpc {
    static const LinearSolverEnum solverType = Gcpc;
    static const std::set< Renumbering > setOfAllowedRenumbering;
    static const bool isHPCCompliant = false;
    static const Renumbering defaultRenumbering = Sans;
};

/**
 * @struct RenumberingChecker
 * @brief Struct statiquepermetant de verifier si un renumeroteur est autorise
         pour un solveur donne
 * @author Nicolas Sellenet
 */
template < class Wrapping > struct RenumberingChecker {
    static bool isHPCCompliant() { return Wrapping::isHPCCompliant; }
};

/** @typedef Definition du verificateur de renumeroteur pour MultFront */
typedef RenumberingChecker< WrapMultFront > MultFrontRenumberingChecker;
/** @typedef Definition du verificateur de renumeroteur pour Ldlt */
typedef RenumberingChecker< WrapLdlt > LdltRenumberingChecker;
/** @typedef Definition du verificateur de renumeroteur pour Mumps */
typedef RenumberingChecker< WrapMumps > MumpsRenumberingChecker;
/** @typedef Definition du verificateur de renumeroteur pour Petsc */
typedef RenumberingChecker< WrapPetsc > PetscRenumberingChecker;
/** @typedef Definition du verificateur de renumeroteur pour Gcpc */
typedef RenumberingChecker< WrapGcpc > GcpcRenumberingChecker;

/**
 * @class BaseLinearSolver
 * @brief Cette classe permet de definir un solveur lineaire
 * @author Nicolas Sellenet
 * @todo verifier que tous les mots-clés sont modifiables par des set
 */
class BaseLinearSolver : public DataStructure {
  protected:
    /** @brief Type du solveur lineaire */
    LinearSolverEnum _linearSolver;
    /** @brief Le solveur est-il vide ? */
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
     * @typedef BaseLinearSolverPtr
     * @brief Pointeur intelligent vers un BaseLinearSolver
     */
    typedef boost::shared_ptr< BaseLinearSolver > BaseLinearSolverPtr;

    /**
     * @brief Constructeur
     * @param currentBaseLinearSolver Type de solveur
     * @param currentRenumber Type de renumeroteur
     * @todo recuperer le code retour de isAllowedRenumberingForSolver
     */
    BaseLinearSolver( const LinearSolverEnum currentBaseLinearSolver = MultFront )
        : BaseLinearSolver( ResultNaming::getNewResultName(), currentBaseLinearSolver ){};

    /**
     * @brief Constructeur
     * @param name Name of the DataStructure
     * @param currentBaseLinearSolver Type de solveur
     * @param currentRenumber Type de renumeroteur
     * @todo recuperer le code retour de isAllowedRenumberingForSolver
     */
    BaseLinearSolver( const std::string name,
                      const LinearSolverEnum currentBaseLinearSolver = MultFront );

    /**
     * @brief Destructor
     */
    ~BaseLinearSolver() { Py_XDECREF( _keywords ); };

    // /** @brief Returns a ListSyntaxMapContainer object "listsyntax",
    //     ready to be inserted  in a CommandSyntax object with the key SOLVEUR
    // */
    // ListSyntaxMapContainer buildListSyntax();

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
     * @brief Methode permettant de savoir si le HPC est autorise
     * @return true si le découpage de domain est autorisé
     */
    virtual bool isHPCCompliant() { return false; };

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
 * @typedef BaseLinearSolverPtr
 * @brief Pointeur intelligent vers un BaseLinearSolver
 */
typedef boost::shared_ptr< BaseLinearSolver > BaseLinearSolverPtr;

/**
 * @class LinearSolver
 * @brief Cette classe permet de definir un solveur lineaire
 * @author Nicolas Sellenet
 * @todo verifier que tous les mots-clés sont modifiables par des set
 */
template < typename linSolvWrap > class LinearSolver : public BaseLinearSolver {
  public:
    /**
     * @typedef LinearSolverPtr
     * @brief Pointeur intelligent vers un LinearSolver
     */
    typedef boost::shared_ptr< LinearSolver< linSolvWrap > > LinearSolverPtr;

    /**
     * @brief Constructeur
     */
    LinearSolver() : LinearSolver( ResultNaming::getNewResultName() ){};

    LinearSolver( const std::string name ) : BaseLinearSolver( name, linSolvWrap::solverType ){};

    bool isHPCCompliant() { return RenumberingChecker< linSolvWrap >::isHPCCompliant(); };
};

typedef LinearSolver< WrapMultFront > MultFrontSolver;
typedef LinearSolver< WrapLdlt > LdltSolver;
typedef LinearSolver< WrapMumps > MumpsSolver;
typedef LinearSolver< WrapPetsc > PetscSolver;
typedef LinearSolver< WrapGcpc > GcpcSolver;

/** @brief Enveloppe d'un pointeur intelligent vers un BaseLinearSolver< MultFront > */
typedef boost::shared_ptr< MultFrontSolver > MultFrontSolverPtr;

/** @brief Enveloppe d'un pointeur intelligent vers un BaseLinearSolver< Ldlt > */
typedef boost::shared_ptr< LdltSolver > LdltSolverPtr;

/** @brief Enveloppe d'un pointeur intelligent vers un BaseLinearSolver< Mumps > */
typedef boost::shared_ptr< MumpsSolver > MumpsSolverPtr;

/** @brief Enveloppe d'un pointeur intelligent vers un BaseLinearSolver< Petsc > */
typedef boost::shared_ptr< PetscSolver > PetscSolverPtr;

/** @brief Enveloppe d'un pointeur intelligent vers un BaseLinearSolver< Gcpc > */
typedef boost::shared_ptr< GcpcSolver > GcpcSolverPtr;

#endif /* LINEARSOLVER_H_ */
