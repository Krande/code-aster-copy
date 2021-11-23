#ifndef DISCRETEPROBLEM_H_
#define DISCRETEPROBLEM_H_

/**
 * @file DiscreteComputation.h
 * @brief Fichier entete de la classe DiscreteComputation
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

#include "Behaviours/BehaviourProperty.h"
#include "LinearAlgebra/ElementaryMatrix.h"
#include "LinearAlgebra/ElementaryVector.h"
#include "Numbering/DOFNumbering.h"
#include "Studies/PhysicalProblem.h"
#include <vector>

/**
 * @class DiscreteComputation
 * @brief Cette classe permet de definir une étude au sens Aster
 * @author Nicolas Sellenet
 */
class DiscreteComputation {
  private:
    /** @brief Etude definie par l'utilisateur */
    PhysicalProblemPtr _study;

    /**
     * @brief Production d'un CommandSyntax pour CALC_MATR_ELEM
     */
    SyntaxMapContainer computeMatrixSyntax( const std::string &optionName );

    /**
     * @brief Calcul des matrices elementaires pour une option quelconque
     */
    ElementaryMatrixDisplacementRealPtr computeMechanicalMatrix( const std::string &optionName );

  public:
    /**
     * @typedef DiscreteComputationPtr
     * @brief Pointeur intelligent vers un DiscreteComputation
     */
    typedef boost::shared_ptr< DiscreteComputation > DiscreteComputationPtr;

    DiscreteComputation( void ) = delete;

    /**
     * @brief Constructeur
     * @param PhysicalProblemPtr Etude utilisateur
     */
    DiscreteComputation( const PhysicalProblemPtr &currentStudy ) : _study( currentStudy ){};

    /**
     * @brief Desctructeur
     */
    ~DiscreteComputation(){};

    /**
     * @brief Fonction permettant de calculer les vecteurs des
              chargements de Dirichlet B*U_imp
     * @param time Instant de calcul
     * @return Vecteur assemblé de chargement
     */
    FieldOnNodesRealPtr imposedDisplacement( ASTERDOUBLE time = 0. );

    /**
     * @brief Fonction permettant de calculer le vecteur pour les
              réactions de Dirichlet B^T * \lambda
     * @param time Instant de calcul
     * @return Vecteur de réaction assemblé
     */
    FieldOnNodesRealPtr dualReaction( FieldOnNodesRealPtr lagr_curr );

    /**
     * @brief Fonction permettant de calculer le vecteur pour les
              déplacements de Dirichlet B * U
     * @param time Instant de calcul
     * @return Vecteur de déplacements de Dirichlet assemblé
     */
    FieldOnNodesRealPtr dualDisplacement( FieldOnNodesRealPtr disp_curr,
                                          ASTERDOUBLE scaling = 1.0 );

    /**
     * @brief Fonction permettant de calculer les vecteurs  pour les
              chargements de Neumann
     * @param time Instants de calcul (vecteur de longueur 3 : instant courant, deltat, paramètre
     theta
     * @return Vecteur des chargement de Neumann assemblé
     */
    FieldOnNodesRealPtr neumann( const VectorReal time, ExternalStateVariablesBuilderPtr );

    /**
     * @brief Fonction permettant de calculer les matrices élémentaires de rigidité
     * @param time Instant de calcul
     * @return Vecteur élémentaire contenant la rigidité mécanique
     */
    ElementaryMatrixDisplacementRealPtr computeElementaryStiffnessMatrix( ASTERDOUBLE time = 0. );

    /**
     * @brief Fonction permettant de calculer les matrices élémentaires pour la matrice tangente
     * utilisée pour l'étape de prédiction de la méthode de Newton
     * @param time Instant de calcul
     * @return Matrice élémentaire contenant la rigidité mécanique
     */
    ElementaryMatrixDisplacementRealPtr computeElementaryTangentMatrix( ASTERDOUBLE time = 0. );

    ElementaryMatrixDisplacementRealPtr computeElementaryJacobianMatrix( ASTERDOUBLE time = 0. );

    /**
     * @brief Construction d'un vecteur de chargement cinématique
     * @return Booleen indiquant que tout s'est bien passe
     */
    FieldOnNodesRealPtr dirichletBC( const ASTERDOUBLE &time ) const;

    /**
     * @brief Construction d'un vecteur de chargement cinématique
     * @return Booleen indiquant que tout s'est bien passe
     */
    FieldOnNodesRealPtr incrementalDirichletBC( const ASTERDOUBLE &time,
                                                const FieldOnNodesRealPtr disp_curr ) const;

    /**
     * @brief Calcul des matrices elementaires pour l'option AMOR_MECA
     */
    ElementaryMatrixDisplacementRealPtr
    computeMechanicalDampingMatrix( const ElementaryMatrixDisplacementRealPtr &rigidity,
                                    const ElementaryMatrixDisplacementRealPtr &mass );

    /**
     * @brief Calcul des matrices elementaires pour l'option RIGI_MECA
     */
    ElementaryMatrixDisplacementRealPtr computeMechanicalStiffnessMatrix();

    /**
     * @brief Calcul des matrices elementaires pour l'option MASS_MECA
     */
    ElementaryMatrixDisplacementRealPtr computeMechanicalMassMatrix();

    /**
     * @brief Récupération de l'étude
     * @return Numérotation du problème discret
     */
    PhysicalProblemPtr getPhysicalProblem() const { return _study; };
};

/**
 * @typedef DiscreteComputationPtr
 * @brief Pointeur intelligent vers un DiscreteComputation
 */
typedef boost::shared_ptr< DiscreteComputation > DiscreteComputationPtr;

#endif /* DISCRETEPROBLEM_H_ */
