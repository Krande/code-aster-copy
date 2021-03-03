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

#include <assert.h>
#include <stdexcept>
#include <string>
#include "aster_fort_calcul.h"
#include "Behaviours/BehaviourProperty.h"

/**
 * @class BehaviourPropertyClass
 * @brief Class to define behaviour
 */

/** @brief Create objects (maps) */
void BehaviourPropertyClass :: createObjects( )
{
    _COMPOR = ConstantFieldOnCellsChar16Ptr( 
              new ConstantFieldOnCellsChar16Class( _baseName + ".COMPOR", _mesh ));

    _MULCOM = ConstantFieldOnCellsChar16Ptr( 
              new ConstantFieldOnCellsChar16Class( _baseName + ".MULCOM", _mesh ));

    _CARCRI = ConstantFieldOnCellsRealPtr( 
              new ConstantFieldOnCellsRealClass( _baseName + ".CARCRI", _mesh ));
};

/** @brief Constructor */
BehaviourPropertyClass :: BehaviourPropertyClass(
  ModelPtr         model,
  MaterialFieldPtr materialField ) :
    _initialState( 0 ), 
    _Implex( 0 ),
    _verbosity( 0 ),
    _model( model ),
    _materialField( materialField )
{
    _mesh = _model -> getMesh( ); 
    _baseName = TemporaryDataStructureNaming::getNewTemporaryName( );
    createObjects( );
};

/** @brief Build objects (maps) */
void BehaviourPropertyClass :: buildObjects( )
{
    std::string modelName = getModel( ) -> getName( );
    modelName.resize( 8, ' ' );

    std::string materialFieldName = getMaterialField( ) -> getName( );
    materialFieldName.resize( 8, ' ' );

    std::string mapComporName =  _COMPOR -> getName( );
    mapComporName.resize( 19, ' ' );

    ASTERINTEGER initialState = getInitialState( );
    ASTERINTEGER implex = getImplex( );
    ASTERINTEGER verbosity = getVerbosity( );

    CALLO_NMDOCC_WRAP( modelName, materialFieldName, mapComporName,
                       &initialState, &implex, &verbosity );

    CALLO_NMDOCR_WRAP( getModel( ) -> getName( ), _CARCRI -> getName( ), &implex );

    CALLO_NMDOCM( getModel( ) -> getName( ), _MULCOM -> getName( ) );
};

