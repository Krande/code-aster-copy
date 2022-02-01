/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#include "LinearAlgebra/ElementaryVector.h"

#include "aster_fort_calcul.h"
#include "astercxx.h"

#include "Supervis/CommandSyntax.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Tools.h"

#include <stdexcept>

FieldOnNodesRealPtr ElementaryVector::assemble( const BaseDOFNumberingPtr dofNume ) const {
    if ( _isEmpty )
        raiseAsterError( "The ElementaryVector is empty" );

    if ( ( !dofNume ) || dofNume->isEmpty() )
        raiseAsterError( "Numerotation is empty" );

    FieldOnNodesRealPtr field = boost::make_shared< FieldOnNodesReal >( dofNume );

    VectorString vect_elem( 1, getName() );

    char *tabNames = vectorStringAsFStrArray( vect_elem, 19 );

    ASTERDOUBLE list_coef = 1.0;
    ASTERINTEGER typscal = 1;
    ASTERINTEGER nbElem = 1;

    std::string base( "G" );
    std::string blanc( "        " );
    std::string zero( "ZERO" );

    CALL_ASSVEC( base.c_str(), field->getName().c_str(), &nbElem, tabNames, &list_coef,
                 dofNume->getName().c_str(), blanc.c_str(), zero.c_str(), &typscal );

    FreeStr( tabNames );

    field->build();

    return field;
};

FieldOnNodesRealPtr ElementaryVector::assembleWithLoadFunctions( const BaseDOFNumberingPtr &dofNume,
                                                                 const ASTERDOUBLE &time ) {
    if ( _isEmpty )
        raiseAsterError( "The ElementaryVector is empty" );

    if ( ( !dofNume ) || dofNume->isEmpty() )
        raiseAsterError( "Numerotation is empty" );

    FieldOnNodesRealPtr field = boost::make_shared< FieldOnNodesReal >( dofNume );
    std::string name( " " );
    name.resize( 24, ' ' );

    _listOfLoads->build();
    if ( !_corichRept->exists() ) {
        _listOfElementaryTerms->updateValuePointer();
        auto size = _listOfElementaryTerms->size();
        for ( ASTERINTEGER i = 1; i <= size; ++i ) {
            std::string vectElem( ( *_listOfElementaryTerms )[i - 1].toString() );
            vectElem.resize( 24, ' ' );
            CALLO_CORICHWRITE( vectElem, &i );
        }
    }

    std::string typres( "R" );
    CALLO_ASASVE( getName(), dofNume->getName(), typres, name );

    std::string detr( "D" );
    std::string fomult( " " );
    const JeveuxVectorChar24 lOF = _listOfLoads->getListOfFunctions();
    if ( !lOF.isEmpty() )
        fomult = lOF->getName();
    std::string param( "INST" );

    JeveuxVectorChar24 vectTmp2( name );
    vectTmp2->updateValuePointer();
    std::string name2( ( *vectTmp2 )[0].toString(), 0, 19 );
    FieldOnNodesRealPtr vectTmp3( new FieldOnNodesReal( name2 ) );
    field->allocateFrom( *vectTmp3 );
    std::string base = "G";

    CALLO_ASCOVA( detr, name, fomult, param, &time, typres, field->getName(), base );

    field->build();

    return field;
};

bool ElementaryVector::build() {
    _listOfElementaryTerms->updateValuePointer();
    _realVector.clear();
    auto size = _listOfElementaryTerms->size();
    _realVector.reserve( size );
    for ( int pos = 0; pos < size; ++pos ) {
        const std::string name = ( *_listOfElementaryTerms )[pos].toString();
        if ( trim( name ) != "" ) {
            ElementaryTermRealPtr toPush =
                boost::make_shared< ElementaryTerm< ASTERDOUBLE > >( name );
            _realVector.push_back( toPush );
        }
    }
    return true;
};
