/**
 * @file BehaviourProperty.cxx
 * @brief Implementation for class BehaviourProperty
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

#include "Behaviours/BehaviourProperty.h"
#include "aster_fort_calcul.h"
#include <assert.h>
#include <stdexcept>
#include <string>

/**
 * @class BehaviourProperty
 * @brief Class to define behaviour
 */

/** @brief Create objects (maps) */
void BehaviourProperty::createObjects() {
    _COMPOR = boost::make_shared< ConstantFieldOnCellsChar16 >( getName() + ".COMPOR    ", _mesh );

    _MULCOM = boost::make_shared< ConstantFieldOnCellsChar16 >( getName() + ".MULCOM", _mesh );

    _CARCRI = boost::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARCRI", _mesh );
};

/** @brief Constructor */
BehaviourProperty::BehaviourProperty()
    : DataStructure( ResultNaming::getNewResultName(), 8, "COMPOR" ), _initialState( false ),
      _implex( false ), _verbosity( false ), _model( nullptr ), _materialField( nullptr ),
      _mesh( nullptr ), _CARCRI( nullptr ), _MULCOM( nullptr ), _COMPOR( nullptr ){};

/** @brief Constructor */
BehaviourProperty::BehaviourProperty( ModelPtr model, MaterialFieldPtr materialField )
    : BehaviourProperty() {
    _model = model;
    _mesh = model->getMesh();
    _materialField = materialField;
};

/** @brief Build objects (maps) */
bool BehaviourProperty::build() {
    createObjects();

    std::string modelName = getModel()->getName();
    modelName.resize( 8, ' ' );

    std::string materialFieldName = getMaterialField()->getName();
    materialFieldName.resize( 8, ' ' );

    std::string comporName = _COMPOR->getName();
    comporName.resize( 19, ' ' );

    std::string base( "G" );

    CALLO_NMDOCC( modelName, materialFieldName, (ASTERLOGICAL *)&_initialState,
                  (ASTERLOGICAL *)&_implex, comporName, base, (ASTERLOGICAL *)&_verbosity );

    CALLO_NMDOCR( getModel()->getName(), _CARCRI->getName(), (ASTERLOGICAL *)&_implex, base );

    CALLO_NMDOCM( getModel()->getName(), _MULCOM->getName(), base );

    AS_ASSERT( _COMPOR->updateValuePointers() );
    AS_ASSERT( _MULCOM->updateValuePointers() );
    AS_ASSERT( _CARCRI->updateValuePointers() );

    return true;
};
