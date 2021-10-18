#ifndef PHYSICALPROBLEM_H_
#define PHYSICALPROBLEM_H_

/**
 * @file PhysicalProblem.h
 * @brief Fichier entete de la classe PhysicalProblem
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

#include "astercxx.h"

#include "Behaviours/BehaviourProperty.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "Loads/ListOfLoads.h"
#include "Materials/BaseExternalStateVariables.h"
#include "Materials/CodedMaterial.h"
#include "Materials/ExternalStateVariablesBuilder.h"
#include "Materials/MaterialField.h"
#include "Modeling/Model.h"
#include "Numbering/DOFNumbering.h"

/**
 * @class PhysicalProblem
 * @brief Cette classe permet de definir une étude au sens Aster
 * @author Nicolas Sellenet
 */
class PhysicalProblem {
  private:
    /** @brief Mesh */
    BaseMeshPtr _mesh;
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Materiau affecté */
    MaterialFieldPtr _materialField;
    /** @brief Liste des chargements */
    ListOfLoadsPtr _listOfLoads;
    /** @brief Liste des chargements */
    ElementaryCharacteristicsPtr _elemChara;
    /** @brief coded material */
    CodedMaterialPtr _codedMater;
    /** @brief Input variables */
    ExternalStateVariablesBuilderPtr _varCom;
    /** @brief Behaviour properties */
    BehaviourPropertyPtr _behavProp;
    /** @brief Numbering */
    BaseDOFNumberingPtr _dofNume;

  public:
    // No default constructor
    PhysicalProblem( void ) = delete;

    /**
     * @brief Constructeur
     * @param ModelPtr Modèle de l'étude
     * @param MaterialFieldPtr Matériau de l'étude
     */
    PhysicalProblem( const ModelPtr curModel, const MaterialFieldPtr curMat,
                     const ElementaryCharacteristicsPtr cara = nullptr );

    ~PhysicalProblem(){};

    /**
     * @brief Add a load (mechanical or dirichlet) with function, formula...
     * @param Args... template list of arguments (load and function or formula)
     */
    template < typename... Args > void addLoad( const Args &...a ) {
        _listOfLoads->addLoad( a... );
    };

    /**
     * @brief Obtenir le matériau affecté
     */
    MaterialFieldPtr getMaterialField() const { return _materialField; };

    /**
     * @brief Obtenir le modèle de l'étude
     */
    ModelPtr getModel() const { return _model; };

    /**
     * @brief Obtenir le maillage de l'étude
     */
    BaseMeshPtr getMesh() const { return _mesh; };
    /**
     * @brief Get elementary characteristics
     */
    ElementaryCharacteristicsPtr getElementaryCharacteristics() const { return _elemChara; };

    /**
     * @brief Get the build coded material
     */
    CodedMaterialPtr getCodedMaterial() const { return _codedMater; };

    /**
     * @brief Renvoit la liste de chargements
     */
    ListOfLoadsPtr getListOfLoads() const { return _listOfLoads; };

    /**
     * @brief Renvoit la carte COMPOR
     */
    BehaviourPropertyPtr getBehaviourProperty() const { return _behavProp; };

    /**
     * @brief Renvoit la carte COMPOR
     */
    BaseDOFNumberingPtr getDOFNumbering() const { return _dofNume; };

    /**
     * @brief Renvoit la carte COMPOR
     */
    void setDOFNumbering( const BaseDOFNumberingPtr dofNume );

    bool computeDOFNumbering();

    /**
     * @brief Create ConstantFieldOnCell for behaviours
     */
    void computeBehaviourProperty( PyObject *keywords, const std::string &initialState,
                                   const std::string &implex, const ASTERINTEGER verbosity );

    void computeBehaviourProperty( PyObject *keywords ) {
        computeBehaviourProperty( keywords, "NON", "NON", 1 );
    };

    /**
     * @brief Construction de la liste de chargements
     * @return true si tout s'est bien passé
     */
    bool computeListOfLoads() { return _listOfLoads->build( _model ); };
};

/**
 * @typedef PhysicalProblemPtr
 * @brief Pointeur intelligent vers un PhysicalProblem
 */
typedef boost::shared_ptr< PhysicalProblem > PhysicalProblemPtr;

#endif /* PHYSICALPROBLEM_H_ */
