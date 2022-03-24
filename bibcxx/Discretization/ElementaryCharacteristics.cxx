/**
 * @file ElementaryCharacteristics.cxx
 * @brief Implementation de ElementaryCharacteristics
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

#include "Discretization/ElementaryCharacteristics.h"

#include "astercxx.h"

#include "DataFields/FieldBuilder.h"

ElementaryCharacteristics::ElementaryCharacteristics( const std::string name,
                                                      const ModelPtr &model )
    : DataStructure( name, 8, "CARA_ELEM" ),
      _model( model ),
      _mesh( model->getMesh() ),
      _CARORIEN( new ConstantFieldOnCellsReal( getName() + ".CARORIEN", _mesh ) ),
      _CARDISCK( new ConstantFieldOnCellsReal( getName() + ".CARDISCK", _mesh ) ),
      _CARDISCM( new ConstantFieldOnCellsReal( getName() + ".CARDISCM", _mesh ) ),
      _CARDISCA( new ConstantFieldOnCellsReal( getName() + ".CARDISCA", _mesh ) ),
      _CARGEOPO( new ConstantFieldOnCellsReal( getName() + ".CARGEOPO", _mesh ) ),
      _CARGENPO( new ConstantFieldOnCellsReal( getName() + ".CARGENPO", _mesh ) ),
      _CARCOQUE( new ConstantFieldOnCellsReal( getName() + ".CARCOQUE", _mesh ) ),
      _CARARCPO( new ConstantFieldOnCellsReal( getName() + ".CARARCPO", _mesh ) ),
      _CARCABLE( new ConstantFieldOnCellsReal( getName() + ".CARCABLE", _mesh ) ),
      _CARGENBA( new ConstantFieldOnCellsReal( getName() + ".CARGENBA", _mesh ) ),
      _CARMASSI( new ConstantFieldOnCellsReal( getName() + ".CARMASSI", _mesh ) ),
      _CARPOUFL( new ConstantFieldOnCellsReal( getName() + ".CARPOUFL", _mesh ) ),
      _CANBSP( new FieldOnCellsLong( getName() + ".CANBSP" ) ),
      _CAFIBR( new FieldOnCellsReal( getName() + ".CAFIBR" ) ),
      _CARDINFO( new ConstantFieldOnCellsReal( getName() + ".CARDINFO", _mesh ) ),
      _model_name( JeveuxVectorChar8( getName() + ".MODELE" ) ),
      _lineic( new ConstantFieldOnCellsChar8( getName() + ".CVENTCXF", _mesh ) ),
      _infos( new ConstantFieldOnCellsReal( getName() + ".CARDINFO", _mesh ) ),
      _isEmpty( true ){};

ModelPtr ElementaryCharacteristics::getModel() const {
    if ( !_model )
        throw std::runtime_error( "Model is empty" );
    return _model;
};

BaseMeshPtr ElementaryCharacteristics::getMesh() const {
    AS_ASSERT( _mesh );

    return _mesh;
};
