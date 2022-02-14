/**
 * @file ElementaryCompute.cxx
 * @brief Implementation of class ElementaryCompute
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

#include "Discretization/ElementaryCompute.h"

void ElementaryCompute::createDescriptor( const ModelPtr &currModel,
                                          const MaterialFieldPtr &currMaterialField,
                                          const ElementaryCharacteristicsPtr &currElemChara ) {
    _rerr->reserve( 5 );
    _rerr->push_back( currModel->getName() );
    _rerr->push_back( _option );
    if ( currModel->nbSuperElement() == 0 ) {
        _rerr->push_back( "NON_SOUS_STRUC" );
    } else {
        _rerr->push_back( "OUI_SOUS_STRUC" );
    }
    _rerr->push_back( currMaterialField->getName() );
    if ( currElemChara != nullptr ) {
        _rerr->push_back( currElemChara->getName() );
    }
};

void ElementaryCompute::createListOfElementaryTerms() { _relr->reserve( 10 ); }
