/**
 * @file BehaviourProperty.cxx
 * @brief Implementation for class BehaviourProperty
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

#include "Behaviours/BehaviourProperty.h"

#include "aster_fort_calcul.h"

#include <stdexcept>
#include <string>

#include <assert.h>

/**
 * @class BehaviourProperty
 * @brief Class to define behaviour
 */

/** @brief Create objects (maps) */
void BehaviourProperty::createObjects() {
    _COMPOR = std::make_shared< ConstantFieldOnCellsChar16 >( getName() + ".COMPOR", _mesh );

    _MULCOM = std::make_shared< ConstantFieldOnCellsChar16 >( getName() + ".MULCOM", _mesh );

    _CARCRI = std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARCRI", _mesh );
};

/** @brief Constructor */
BehaviourProperty::BehaviourProperty( const std::string name )
    : DataStructure( name, 8, "COMPOR" ),
      _initialState( false ),
      _verbosity( false ),
      _model( nullptr ),
      _materialField( nullptr ),
      _mesh( nullptr ),
      _CARCRI( nullptr ),
      _MULCOM( nullptr ),
      _COMPOR( nullptr ),
      _annealing( false ) {};

/** @brief Constructor */
BehaviourProperty::BehaviourProperty() : BehaviourProperty( ResultNaming::getNewResultName() ) {};

/** @brief Constructor */
BehaviourProperty::BehaviourProperty( const std::string name, ModelPtr model,
                                      MaterialFieldPtr materialField )
    : BehaviourProperty( name ) {
    _model = model;
    _mesh = model->getMesh();
    _materialField = materialField;
};

/** @brief Constructor */
BehaviourProperty::BehaviourProperty( ModelPtr model, MaterialFieldPtr materialField )
    : BehaviourProperty( ResultNaming::getNewResultName(), model, materialField ) {};

/** @brief Build objects (maps) */
bool BehaviourProperty::build() {
    createObjects();

    std::string modelName = getModel()->getName();
    modelName.resize( 8, ' ' );

    std::string comporName = _COMPOR->getName();
    comporName.resize( 19, ' ' );

    std::string base( "G" );

    if ( getModel()->isMechanical() ) {

        std::string materialFieldName = getMaterialField()->getName();
        materialFieldName.resize( 8, ' ' );

        CALLO_NMDOCC( modelName, materialFieldName, (ASTERLOGICAL *)&_initialState, comporName,
                      base, (ASTERLOGICAL *)&_verbosity );

        CALLO_NMDOCR( getModel()->getName(), _CARCRI->getName(), base );

        CALLO_NMDOCM( getModel()->getName(), _MULCOM->getName(), base );

        _MULCOM->updateValuePointers();
        _CARCRI->updateValuePointers();
        detectFunctionnalities();
    } else if ( getModel()->isThermal() ) {
        CALLO_NXDOCC( modelName, comporName, base );
    } else {
        AS_ABORT( "Not implemented" );
    }

    _COMPOR->updateValuePointers();

    return true;
};

bool BehaviourProperty::hasBehaviour( const std::string &behaviour ) const {

    if ( _COMPOR && _COMPOR->exists() ) {
        auto values = _COMPOR->getValues();
        for ( auto &zone : values ) {
            auto cmps = zone.getComponents();
            auto val = zone.getValues();
            for ( int i = 0; i < cmps.size(); i++ ) {
                if ( strip( cmps[i] ) == "RELCOM" ) {
                    if ( val[i].toString() == behaviour )
                        return true;
                }
            }
        }
    }
    return false;
};

void BehaviourProperty::detectFunctionnalities() {
    // Detect annealing
    std::string feature( "Annealing" );
    feature.resize( 16, ' ' );
    std::string caraElem( " " );
    caraElem.resize( 8, ' ' );
    bool flag;
    CALLO_HASBEHAVIOURFEATURE( this->getModel()->getName(), caraElem,
                               this->getBehaviourField()->getName(), feature,
                               (ASTERLOGICAL *)&flag );
    _annealing = false;
    if ( flag ) {
        _annealing = true;
    }
};
