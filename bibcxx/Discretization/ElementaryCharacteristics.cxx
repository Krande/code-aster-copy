/**
 * @file ElementaryCharacteristics.cxx
 * @brief Implementation de ElementaryCharacteristics
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

#include "Discretization/ElementaryCharacteristics.h"

#include "DataFields/FieldBuilder.h"

ElementaryCharacteristics::ElementaryCharacteristics( const std::string name,
                                                      const ModelPtr &model )
    : DataStructure( name, 8, "CARA_ELEM" ),
      _model( model ),
      _mesh( model->getMesh() ),
      _CARORIEN( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARORIEN", _mesh ) ),
      _CARDISCK( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARDISCK", _mesh ) ),
      _CARDISCM( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARDISCM", _mesh ) ),
      _CARDISCA( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARDISCA", _mesh ) ),
      _CARGEOPO( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARGEOPO", _mesh ) ),
      _CARGENPO( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARGENPO", _mesh ) ),
      _CARCOQUE( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARCOQUE", _mesh ) ),
      _CARARCPO( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARARCPO", _mesh ) ),
      _CARCABLE( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARCABLE", _mesh ) ),
      _CARGENBA( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARGENBA", _mesh ) ),
      _CARMASSI( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARMASSI", _mesh ) ),
      _CARPOUFL( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARPOUFL", _mesh ) ),
      _CANBSP( std::make_shared< FieldOnCellsLong >( getName() + ".CANBSP" ) ),
      _CAFIBR( std::make_shared< FieldOnCellsReal >( getName() + ".CAFIBR" ) ),
      _CARDINFO( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARDINFO", _mesh ) ),
      _model_name( JeveuxVectorChar8( getName() + ".MODELE" ) ),
      _lineic( std::make_shared< ConstantFieldOnCellsChar8 >( getName() + ".CVENTCXF", _mesh ) ),
      _infos( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CARDINFO", _mesh ) ),
      _isEmpty( true ) {};

ModelPtr ElementaryCharacteristics::getModel() const { return _model; };

BaseMeshPtr ElementaryCharacteristics::getMesh() const { return _mesh; };
