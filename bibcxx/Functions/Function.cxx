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

#include "Functions/Function.h"

#include "Supervis/ResultNaming.h"

#include <stdexcept>
#include <string>
#include <vector>

FunctionPtr emptyRealFunction( new Function( "" ) );

BaseFunction::BaseFunction( const std::string name, const std::string type,
                            const std::string type2 )
    : GenericFunction( name, type, type2 ), _value( JeveuxVectorReal( getName() + ".VALE" ) ) {}

BaseFunction::BaseFunction( const std::string type, const std::string type2 )
    : BaseFunction::BaseFunction( ResultNaming::getNewResultName(), type, type2 ) {}

void BaseFunction::allocate( ASTERINTEGER size ) {
    if ( _property.exists() )
        _property->deallocate();
    propertyAllocate();

    if ( _value.exists() )
        _value->deallocate();
    _value->allocate( 2 * size );
}

void BaseFunction::deallocate() {
    _property->deallocate();
    _value->deallocate();
}

std::string BaseFunction::getResultName() {
    if ( !_property.exists() )
        return "";
    _property->updateValuePointer();
    return ( *_property )[3].toString();
}

void BaseFunction::setParameterName( const std::string name ) {
    if ( !_property.exists() )
        propertyAllocate();
    ( *_property )[2] = name.substr( 0, 8 ).c_str();
}

void BaseFunction::setResultName( const std::string name ) {
    if ( !_property.exists() )
        propertyAllocate();
    ( *_property )[3] = name.substr( 0, 8 ).c_str();
}

void BaseFunction::setValues( const VectorReal &absc, const VectorReal &ordo ) {
    if ( absc.size() != ordo.size() )
        throw std::runtime_error( "Function: length of abscissa and ordinates must be equal" );

    // Create Jeveux vector ".VALE"
    const int nbpts = absc.size();
    if ( _value.exists() && _value->size() != 2 * nbpts )
        throw std::runtime_error( "Function: the function size is " + std::to_string( nbpts ) +
                                  ", lists of this size are expected" );

    if ( !_value.exists() )
        _value->allocate( 2 * nbpts );

    // Loop on the points
    VectorReal::const_iterator abscIt = absc.begin();
    VectorReal::const_iterator ordoIt = ordo.begin();
    int idx = 0;
    for ( ; abscIt != absc.end(); ++abscIt, ++ordoIt ) {
        ( *_value )[idx] = *abscIt;
        ( *_value )[nbpts + idx] = *ordoIt;
        ++idx;
    }
}

void BaseFunction::setInterpolation( const std::string type ) {
    std::string interp;
    if ( !_property.exists() )
        propertyAllocate();

    if ( type.length() != 7 )
        throw std::runtime_error( "Function: interpolation must be 7 characters long." );

    interp = type.substr( 0, 3 );
    if ( interp != "LIN" && interp != "LOG" && interp != "NON" )
        throw std::runtime_error( "Function: invalid interpolation for abscissa." );

    interp = type.substr( 4, 3 );
    if ( interp != "LIN" && interp != "LOG" && interp != "NON" )
        throw std::runtime_error( "Function: invalid interpolation for ordinates." );

    ( *_property )[1] = type.c_str();
}

void BaseFunction::setAsConstant() {
    if ( !_property.exists() )
        propertyAllocate();
    _funct_type = "CONSTANT";
    ( *_property )[0] = _funct_type;
}

/* Complex function */
void FunctionComplex::allocate( ASTERINTEGER size ) {
    if ( _property.exists() )
        _property->deallocate();
    propertyAllocate();

    if ( _value.exists() )
        _value->deallocate();
    _value->allocate( 3 * size );
}

void FunctionComplex::setValues( const VectorReal &absc, const VectorReal &ordo ) {
    if ( absc.size() * 2 != ordo.size() )
        throw std::runtime_error(
            "Function: The length of ordinates must be twice that of abscissas." );

    // Create Jeveux vector ".VALE"
    const int nbpts = absc.size();
    if ( _value.exists() && _value->size() != 3 * nbpts )
        throw std::runtime_error( "Function: the function size is " + std::to_string( nbpts ) +
                                  ", abscissas of this size are expected" );

    if ( !_value.exists() )
        _value->allocate( 3 * nbpts );

    // Loop on the points
    VectorReal::const_iterator abscIt = absc.begin();
    VectorReal::const_iterator ordoIt = ordo.begin();
    int idx = 0;
    for ( ; abscIt != absc.end(); ++abscIt, ++ordoIt ) {
        ( *_value )[idx] = *abscIt;
        ( *_value )[nbpts + 2 * idx] = *ordoIt;
        ++ordoIt;
        ( *_value )[nbpts + 2 * idx + 1] = *ordoIt;
        ++idx;
    }
}

void FunctionComplex::setValues( const VectorReal &absc, const VectorComplex &ordo ) {
    throw std::runtime_error( "Not yet implemented!" );
}
