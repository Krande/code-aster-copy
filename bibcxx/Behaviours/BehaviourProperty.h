#ifndef BEHAVIOURPROPERTY_H_
#define BEHAVIOURPROPERTY_H_

/**
 * @file BehaviourProperty.h
 * @brief Header for class BehaviourProperty
 * @author Mickaël Abas
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

/* person_in_charge: mickael.abbas at edf.fr */

#include "astercxx.h"
#include <string>
#include "Modeling/Model.h"
#include "Materials/MaterialField.h"
#include "DataFields/ConstantFieldOnCells.h"
#include "DataStructures/TemporaryDataStructureNaming.h"

// TODO doit être une DataStructure
/**
 * @class BehaviourProperty
 * @brief Class to define behaviour
 */
class BehaviourProperty {
  private:
    /** @brief Mesh */
    BaseMeshPtr _mesh;

    /** @brief Model */
    ModelPtr _model;

    /** @brief Material field */
    MaterialFieldPtr _materialField;

    /** @brief Base name of the maps */
    std::string _baseName;

    /** @brief Flag for initial state */
    bool _initialState;

    /** @brief Flag for IMPLEX algorithm */
    bool _implex;

    /** @brief Verbosity */
    bool _verbosity;

    /** @brief Map '.COMPOR' to define behaviours */
    ConstantFieldOnCellsChar16Ptr _COMPOR;

    /** @brief Mao '.MULCOM' for multi-behaviours (crystals) */
    ConstantFieldOnCellsChar16Ptr _MULCOM;

    /** @brief Map '.CARCRI' for parameters to integrate behaviours */
    ConstantFieldOnCellsRealPtr _CARCRI;

  private:
    /** @brief Create objects (maps) */
    void createObjects( );

  public:
    /** @brief Constructor */
    BehaviourProperty(
      ModelPtr         model,
      MaterialFieldPtr materialField );

    /** @brief Destructor */
    ~BehaviourProperty( ){};

    /** @brief Build object */
    void build( );

    /** @brief Get model */
    ModelPtr getModel( ) const { return _model; }

    /** @brief Get material field */
    MaterialFieldPtr getMaterialField( ) const { return _materialField; }

    /** @brief Set flag for initial state */
    void setInitialState( const bool &value ) { _initialState = value; };

    /** @brief Set flag for Implex */
    void setImplex( const bool &value ) { _implex = value; };

    /** @brief Set verbosity */
    void setVerbosity( const bool &value ) { _verbosity = value; };

};

/** @typedef Smart-pointer to behaviour class */
typedef boost::shared_ptr< BehaviourProperty > BehaviourPropertyPtr;


#endif
