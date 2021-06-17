/**
 * @file BehaviourProperty.cxx
 * @brief Implementation for class BehaviourProperty
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
void BehaviourProperty ::createObjects() {
    _COMPOR = ConstantFieldOnCellsChar16Ptr(
        new ConstantFieldOnCellsChar16( _baseName + ".COMPOR", _mesh ) );

    _MULCOM = ConstantFieldOnCellsChar16Ptr(
        new ConstantFieldOnCellsChar16( _baseName + ".MULCOM", _mesh ) );

    _CARCRI = ConstantFieldOnCellsRealPtr(
        new ConstantFieldOnCellsReal( _baseName + ".CARCRI", _mesh ) );
};

/** @brief Constructor */
BehaviourProperty ::BehaviourProperty( ModelPtr model, MaterialFieldPtr materialField )
    : _initialState( false ), _implex( false ), _verbosity( false ), _model( model ),
      _materialField( materialField ) {
    _mesh = _model->getMesh();
    _baseName = TemporaryDataStructureNaming::getNewTemporaryName();
    createObjects();
};

/** @brief Build objects (maps) */
void BehaviourProperty ::build() {
    std::string modelName = getModel()->getName();
    modelName.resize( 8, ' ' );

    std::string materialFieldName = getMaterialField()->getName();
    materialFieldName.resize( 8, ' ' );

    std::string mapComporName = _COMPOR->getName();
    mapComporName.resize( 19, ' ' );

    CALLO_NMDOCC( modelName, materialFieldName, (ASTERLOGICAL *)&_initialState,
                  (ASTERLOGICAL *)&_implex, mapComporName, (ASTERLOGICAL *)&_verbosity );

    CALLO_NMDOCR( getModel()->getName(), _CARCRI->getName(), (ASTERLOGICAL *)&_implex );

    CALLO_NMDOCM( getModel()->getName(), _MULCOM->getName() );
};
