/**
 * @file Function2D.cxx
 * @brief Implementation de Function2D
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "Functions/Function2D.h"

std::string Function2D::getResultName() {
    if ( !_property.exists() )
        return "";
    _property->updateValuePointer();
    return ( *_property )[3].toString();
}

ASTERINTEGER Function2D::maximumSize() const {
    if ( !_value.exists() )
        return 0;
    _value->build();
    ASTERINTEGER toReturn = 0;
    for ( const auto &obj : *_value ) {
        if ( obj->size() > toReturn )
            toReturn = obj->size();
    }
    return toReturn;
}

ASTERINTEGER Function2D::size() const {
    if ( !_value.exists() )
        return 0;
    _value->build();
    ASTERINTEGER toReturn = 0;
    for ( const auto &obj : *_value ) {
        toReturn += obj->size();
    }
    return toReturn;
}
