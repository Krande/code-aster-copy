/**
 * @file ElementaryCompute.cxx
 * @brief Implementation of class ElementaryCompute
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

void ElementaryCompute::createDescriptor( const ModelPtr &currModel, const std::string &option ) {
    _rerr->allocate( 3 );
    _model = currModel;
    if ( currModel ) {
        ( *_rerr )[0] = currModel->getName();

        if ( currModel->numberOfSuperElement() == 0 ) {
            ( *_rerr )[2] = "NON_SOUS_STRUC";
        } else {
            ( *_rerr )[2] = "OUI_SOUS_STRUC";
        }
    } else {
        ( *_rerr )[2] = "NON_SOUS_STRUC";
    }

    ( *_rerr )[1] = option;
};

void ElementaryCompute::createDescriptor( const ModelPtr &currModel ) {
    if ( _rerr->exists() ) {
        _rerr->updateValuePointer();
    } else {
        _rerr->allocate( 1 );
    }
    _model = currModel;
    if ( _model ) {
        ( *_rerr )[0] = _model->getName();
    }
};

std::string ElementaryCompute::getOption() const {
    if ( _rerr->size() > 1 ) {
        _rerr->updateValuePointer();
        return strip( ( *_rerr )[1].toString() );
    }

    return std::string();
}
