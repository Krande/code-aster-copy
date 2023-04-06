/**
 * @file ResultNaming.cxx
 * @brief Implementation of automatic naming of jeveux objects.
 * @section LICENCE
 * Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
 * This file is part of code_aster.
 *
 * code_aster is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * code_aster is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with code_aster.  If not, see <http://www.gnu.org/licenses/>.

 * person_in_charge: mathieu.courtois@edf.fr
 */

#include "astercxx.h"

#include "Functions/GenericFunction.h"

#include "Supervis/ResultNaming.h"

#include <stdexcept>
#include <string>
#include <vector>

GenericFunction::GenericFunction( const std::string &name, const std::string &type,
                                  const std::string &functType )
    : DataStructure( name, 19, type ),
      _property( JeveuxVectorChar24( getName() + ".PROL" ) ),
      _funct_type( functType ) {}

VectorString GenericFunction::getProperties() const {
    _property->updateValuePointer();
    const auto size = _property->size();
    VectorString prop;
    for ( int i = 0; i < size; ++i ) {
        prop.push_back( ( *_property )[i].rstrip() );
    }
    return prop;
}

void GenericFunction::setExtrapolation( const std::string type ) {
    if ( !_property.exists() )
        propertyAllocate();

    if ( type.length() != 2 )
        throw std::runtime_error( "Function: interpolation must be 2 characters long." );

    std::string auth( "CELI" );
    if ( auth.find( type[0] ) == std::string::npos )
        throw std::runtime_error( "Function: invalid extrapolation for abscissa." );

    if ( auth.find( type[1] ) == std::string::npos )
        throw std::runtime_error( "Function: invalid extrapolation for ordinates." );

    ( *_property )[4] = type.c_str();
}
