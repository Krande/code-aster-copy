#ifndef STUDYDESCRIPTION_H_
#define STUDYDESCRIPTION_H_

/**
 * @file StudyDescription.h
 * @brief Fichier entete de la classe StudyDescription
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

#include "Discretization/ElementaryCharacteristics.h"
#include "Loads/DirichletBC.h"
#include "Loads/ListOfLoads.h"
#include "Loads/MechanicalLoad.h"
#include "Materials/BaseExternalStateVariables.h"
#include "Materials/CodedMaterial.h"
#include "Materials/ExternalStateVariablesBuilder.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Numbering/DOFNumbering.h"

/**
 * @class StudyDescription
 * @brief Cette classe permet de definir une étude au sens Aster
 * @author Nicolas Sellenet
 */
class StudyDescription {
  private:
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

  public:

    // No default constructor
    StudyDescription( void ) = delete;

    /**
     * @brief Constructeur
     * @param ModelPtr Modèle de l'étude
     * @param MaterialFieldPtr Matériau de l'étude
     */
    StudyDescription( const ModelPtr &curModel, const MaterialFieldPtr &curMat,
                              const ElementaryCharacteristicsPtr &cara = nullptr )
        : _model( curModel ), _materialField( curMat ),
          _listOfLoads( boost::make_shared< ListOfLoads >() ), _elemChara( cara ),
          _codedMater( boost::make_shared<  CodedMaterial >( _materialField, _model ) ),
          _varCom( boost::make_shared< ExternalStateVariablesBuilder >( _model, _materialField,
                                                          _elemChara, _codedMater ) ){
        if( _elemChara ){
            if( _model->getName() != _elemChara->getModel()->getName())
                throw std::runtime_error("Inconsistent model");
        }

        if( _model->getMesh()->getName() != _materialField->getMesh()->getName())
            throw std::runtime_error("Inconsistent mesh");
    };

    ~StudyDescription(){};

    /**
     * @brief Add a load (mechanical or dirichlet) with function, formula...
     * @param Args... template list of arguments (load and function or formula)
     */
    template < typename... Args > void addLoad( const Args &... a ) {
        _listOfLoads->addLoad( a... );
    };

    /**
     * @brief Construction de la liste de chargements
     * @return true si tout s'est bien passé
     */
    bool computeListOfLoads() { return _listOfLoads->build(_model); };

    /**
     * @brief Get elementary characteristics
     */
    const ElementaryCharacteristicsPtr &getElementaryCharacteristics() const { return _elemChara; };

    /**
     * @brief Get the build coded material
     */
    const CodedMaterialPtr &getCodedMaterial() const { return _codedMater; };

    /**
     * @brief Obtenir la liste des chargements cinematiques
     */
    const ListDiriBC &getListOfDirichletBCs() const {
        return _listOfLoads->getListOfDirichletBCs();
    };

    /**
     * @brief Renvoit la liste de chargements
     */
    const ListOfLoadsPtr &getListOfLoads() const { return _listOfLoads; };

    /**
     * @brief Obtenir la liste des chargements mecaniques
     */
    const ListMecaLoadReal &getListOfMechanicalLoadsReal() const {
        return _listOfLoads->getListOfMechanicalLoadsReal();
    };

    /**
     * @brief Obtenir la liste des chargements mecaniques
     */
    const ListMecaLoadFunction &getListOfMechanicalLoadsFunction() const {
        return _listOfLoads->getListOfMechanicalLoadsFunction();
    };

#ifdef ASTER_HAVE_MPI
    /**
     * @brief Obtenir la liste des chargements mecaniques
     */
    const ListParaMecaLoadReal &getListOfParallelMechanicalLoadsReal() const {
        return _listOfLoads->getListOfParallelMechanicalLoadsReal();
    };

    /**
     * @brief Obtenir la liste des chargements mecaniques
     */
    const ListParaMecaLoadFunction &getListOfParallelMechanicalLoadsFunction() const {
        return _listOfLoads->getListOfParallelMechanicalLoadsFunction();
    };
#endif /* ASTER_HAVE_MPI */

    /**
     * @brief Function de récupération de la liste des charges thermiques
     * @return _listOfThermalLoads
     */
    const ListTherLoadReal &getListOfThermalLoadsReal() const
    { return _listOfLoads->getListOfThermalLoadsReal(); };

    /**
     * @brief Function de récupération de la liste des charges thermiques
     * @return _listOfThermalLoads
     */
    const ListTherLoadFunction &getListOfThermalLoadsFunction() const
    { return _listOfLoads->getListOfThermalLoadsFunction(); };

    /**
     * @brief Obtenir le matériau affecté
     */
    const MaterialFieldPtr &getMaterialField() const { return _materialField; };

    /**
     * @brief Obtenir le modèle de l'étude
     */
    const ModelPtr &getModel() const { return _model; };

    /**
     * @brief Obtenir le maillage de l'étude
     */
     BaseMeshPtr getMesh() const { return _model->getMesh(); };
};

/**
 * @typedef StudyDescriptionPtr
 * @brief Pointeur intelligent vers un StudyDescription
 */
typedef boost::shared_ptr< StudyDescription > StudyDescriptionPtr;

#endif /* STUDYDESCRIPTION_H_ */
