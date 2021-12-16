#ifndef STATICMECHANICCONTEXT_H_
#define STATICMECHANICCONTEXT_H_

/**
 * @file StaticMechanicalContext.h
 * @brief Fichier entete de la classe StaticMechanicalContext
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

#include "Discretization/DiscreteComputation.h"
#include "Loads/ListOfLoads.h"
#include "Materials/BaseExternalStateVariables.h"
#include "Results/Result.h"
#include "Solvers/LinearSolver.h"

class StaticMechanicalAlgorithm;

/**
 * @class StaticMechanicalContext
 * @brief Context around static mechanical algorithm
 * @author Nicolas Sellenet
 */
class StaticMechanicalContext {
  private:
    /** @brief Problème discret */
    DiscreteComputationPtr _discreteComputation;
    /** @brief Solveur linéaire */
    LinearSolverPtr _linearSolver;
    /** @brief Sd de stockage des résultats */
    ResultPtr _results;
    /** @brief Chargements */
    ListOfLoadsPtr _listOfLoads;
    /** @brief Pas de temps courant */
    ASTERDOUBLE _time;
    /** @brief rank */
    ASTERINTEGER _rank;
    /** @brief Assembly matrix */
    AssemblyMatrixDisplacementRealPtr _aMatrix;
    /** @brief Are elastic properties constant */
    bool _isConst;
    /** @brief Input variables */
    ExternalStateVariablesBuilderPtr _varCom;

  public:
    /** @brief Timer */
    std::map< std::string, ASTERDOUBLE > _timer;

  public:
    /**
     * @brief Constructeur
     * @param DiscreteComputationPtr Problème discret a résoudre par l'algo
     * @param LinearSolverPtr Sovleur linéaire qui sera utilisé
     * @param ResultPtr Résultat pour le stockage des déplacements
     */
    StaticMechanicalContext( const DiscreteComputationPtr &curPb, const LinearSolverPtr linSolv,
                             const ResultPtr container )
        : _discreteComputation( curPb ), _linearSolver( linSolv ),
          _listOfLoads( _discreteComputation->getPhysicalProblem()->getListOfLoads() ),
          _results( container ), _time( 0. ), _rank( 1 ),
          _aMatrix( new AssemblyMatrixDisplacementReal() ),
          _isConst( _discreteComputation->getPhysicalProblem()->getCodedMaterial()->constant() ),
          _varCom( new ExternalStateVariablesBuilder(
              _discreteComputation->getPhysicalProblem()->getModel(),
              _discreteComputation->getPhysicalProblem()->getMaterialField(),
              _discreteComputation->getPhysicalProblem()->getElementaryCharacteristics(),
              _discreteComputation->getPhysicalProblem()->getCodedMaterial() ) ),
          _timer( { { "Matrix", 0.0 },
                    { "Rhs", 0.0 },
                    { "Facto", 0.0 },
                    { "Solve", 0.0 },
                    { "Post", 0.0 } } ){};

    /**
     * @brief Function to set the "position" of the context
     * @param time time value
     * @param rank number of iteration
     */
    void setStep( const ASTERDOUBLE &time, const ASTERINTEGER &rank ) {
        _time = time;
        _rank = rank;
    };

    AssemblyMatrixDisplacementRealPtr getStiffnessMatrix( void ) { return _aMatrix; }

    std::map< std::string, ASTERDOUBLE > getTimer() { return _timer; }

    DiscreteComputationPtr getDiscreteComputation() { return _discreteComputation; }

    ResultPtr getResult() { return _results; }

    friend class StaticMechanicalAlgorithm;
};

#endif /* STATICMECHANICCONTEXT_H_ */
