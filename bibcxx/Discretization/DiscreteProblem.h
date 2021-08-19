#ifndef DISCRETEPROBLEM_H_
#define DISCRETEPROBLEM_H_

/**
 * @file DiscreteProblem.h
 * @brief Fichier entete de la classe DiscreteProblem
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
#include "Studies/StudyDescription.h"
#include <vector>

/**
 * @class DiscreteProblem
 * @brief Cette classe permet de definir une étude au sens Aster
 * @author Nicolas Sellenet
 */
class DiscreteProblem {
  private:
    /** @brief Etude definie par l'utilisateur */
    StudyDescriptionPtr _study;

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
     * @typedef DiscreteProblemPtr
     * @brief Pointeur intelligent vers un DiscreteProblem
     */
    typedef boost::shared_ptr< DiscreteProblem > DiscreteProblemPtr;

    DiscreteProblem( void ) = delete;

    /**
     * @brief Constructeur
     * @param StudyDescriptionPtr Etude utilisateur
     */
    DiscreteProblem( const StudyDescriptionPtr &currentStudy ) : _study( currentStudy ){};

    /**
     * @brief Desctructeur
     */
    ~DiscreteProblem(){};

    /**
     * @brief Calcul des matrices elementaires pour l'option CHAR_MECA
     */
    ElementaryVectorDisplacementRealPtr computeElementaryMechanicalLoadsVector();

    /**
     * @brief Fonction permettant de calculer les vecteurs élémentaires pour les
              chargements de Dirichlet
     * @param time Instant de calcul
     * @return Vecteur élémentaire
     */
    ElementaryVectorDisplacementRealPtr computeElementaryDirichletVector( ASTERDOUBLE time = 0. );

    FieldOnNodesRealPtr
    computeDirichlet( BaseDOFNumberingPtr dofNume,  ASTERDOUBLE time = 0.);

    /**
     * @brief Fonction permettant de calculer les vecteurs élémentaires pour les
              chargements de Dirichlet
     * @param time Instant de calcul
     * @return Vecteur élémentaire
     */
    ElementaryVectorDisplacementRealPtr
    computeElementaryDirichletReactionVector( FieldOnNodesRealPtr lagr_curr );

    FieldOnNodesRealPtr
    computeDirichletReaction( BaseDOFNumberingPtr dofNume,  FieldOnNodesRealPtr lagr_curr);



    ElementaryVectorDisplacementRealPtr
    computeElementaryDualizedDirichletVector( FieldOnNodesRealPtr disp_curr,
                                              ASTERDOUBLE scaling = 1.0 );

    FieldOnNodesRealPtr
    computeDualizedDirichlet( BaseDOFNumberingPtr dofNume,  FieldOnNodesRealPtr disp_curr,
                              ASTERDOUBLE scaling = 1.0);

    /**
     * @brief Fonction permettant de calculer les vecteurs élémentaires pour les
              forces de Laplace
     * @return Vecteur élémentaire
     */
    ElementaryVectorDisplacementRealPtr computeElementaryLaplaceVector();

    /**
     * @brief Fonction permettant de calculer les vecteurs élémentaires pour les
              chargements de Neumann
     * @param time Instants de calcul (vecteur de longueur 3 : instant courant, deltat, paramètre
     theta
     * @return Vecteur élémentaire
     */
    ElementaryVectorDisplacementRealPtr computeElementaryNeumannVector( const VectorReal time,
                                                      ExternalStateVariablesBuilderPtr );

    FieldOnNodesRealPtr
    computeNeumann( BaseDOFNumberingPtr dofNume,  const VectorReal time,
                                                  ExternalStateVariablesBuilderPtr);

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
    FieldOnNodesRealPtr computeDirichletBC( const BaseDOFNumberingPtr &curDOFNum,
                                          const ASTERDOUBLE &time ) const;

    /**
     * @brief Détermination de la numérotation de ddl
     * @return Numérotation du problème discret
     */
    BaseDOFNumberingPtr
    computeDOFNumbering( BaseDOFNumberingPtr dofNum = BaseDOFNumberingPtr( nullptr ) );

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
    StudyDescriptionPtr getStudyDescription() const { return _study; };

    ListOfLoadsPtr getListOfLoads() const { return _study->getListOfLoads(); };

    /**
     * @brief Create ConstantFieldOnCell for behaviours
     */
    BehaviourPropertyPtr createBehaviour( PyObject *keywords,
                          const std::string &initialState = "NON",
                          const std::string &implex = "NON", const int verbosity = 0 );
    // added for Python interface
    BehaviourPropertyPtr createBehaviour( PyObject *keywords )
    { return createBehaviour( keywords, "NON", "NON", 1 ); };
};

/**
 * @typedef DiscreteProblemPtr
 * @brief Pointeur intelligent vers un DiscreteProblem
 */
typedef boost::shared_ptr< DiscreteProblem > DiscreteProblemPtr;

#endif /* DISCRETEPROBLEM_H_ */
