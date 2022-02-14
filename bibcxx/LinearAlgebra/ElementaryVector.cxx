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

FieldOnNodesRealPtr BaseElementaryVector::assemble( const BaseDOFNumberingPtr dofNume ) const {
    if ( _isEmpty )
        raiseAsterError( "The ElementaryVector is empty" );

    if ( ( !dofNume ) || dofNume->isEmpty() )
        raiseAsterError( "Numerotation is empty" );

    // Create field
    FieldOnNodesRealPtr field = boost::make_shared< FieldOnNodesReal >( dofNume );

    // Elementary vector names
    std::string vectElemName = getName();
    VectorString vectElemVect( 1, vectElemName );
    char *tabNames = vectorStringAsFStrArray( vectElemVect, 19 );

    // Assembling
    ASTERDOUBLE list_coef = 1.0;
    ASTERINTEGER typscal = 1;
    ASTERINTEGER nbElem = 1;
    std::string base( "G" );

    CALL_ASSVEC( base.c_str(), field->getName().c_str(), &nbElem, tabNames, &list_coef,
                 dofNume->getName().c_str(), &typscal );

    FreeStr( tabNames );

    return field;
};

FieldOnNodesRealPtr
BaseElementaryVector::assembleWithLoadFunctions( const BaseDOFNumberingPtr &dofNume,
                                                 const ASTERDOUBLE &time ) {
    if ( _isEmpty )
        raiseAsterError( "The ElementaryVector is empty" );

    if ( ( !dofNume ) || dofNume->isEmpty() )
        raiseAsterError( "Numerotation is empty" );

    // Create field
    FieldOnNodesRealPtr field = boost::make_shared< FieldOnNodesReal >( dofNume );

    // Elementary vector names
    std::string vectElemName = getName();

    // Link between vector and load function
    if ( !_corichRept->exists() ) {
        _relr->updateValuePointer();
        auto size = _relr->size();
        for ( ASTERINTEGER i = 1; i <= size; ++i ) {
            std::string vectElem( ( *_relr )[i - 1].toString() );
            vectElem.resize( 24, ' ' );
            CALLO_CORICHWRITE( vectElem, &i );
        }
    }

    // Pre-assembling
    std::string typres( "R" );
    std::string name( " " );
    name.resize( 24, ' ' );
    CALLO_ASASVE( vectElemName, dofNume->getName(), typres, name );

    // Get function for load
    _listOfLoads->build();
    std::string fomult( " " );
    const JeveuxVectorChar24 listOfLoadsFunc = _listOfLoads->getListOfFunctions();
    if ( !listOfLoadsFunc.isEmpty() )
        fomult = listOfLoadsFunc->getName();

    // Final assembling with load function
    JeveuxVectorChar24 vectTmp2( name );
    vectTmp2->updateValuePointer();
    std::string name2( ( *vectTmp2 )[0].toString(), 0, 19 );
    FieldOnNodesRealPtr vectTmp3( new FieldOnNodesReal( name2 ) );
    field->allocateFrom( *vectTmp3 );

    std::string detr( "D" );
    std::string base = "G";
    std::string param( "INST" );
    CALLO_ASCOVA( detr, name, fomult, param, &time, typres, field->getName(), base );

    return field;
};

bool BaseElementaryVector::build() {
    _relr->updateValuePointer();
    _realVector.clear();
    auto size = _relr->size();
    _realVector.reserve( size );
    for ( int pos = 0; pos < size; ++pos ) {
        const std::string name = ( *_relr )[pos].toString();
        if ( trim( name ) != "" ) {
            ElementaryTermRealPtr toPush =
                boost::make_shared< ElementaryTerm< ASTERDOUBLE > >( name );
            _realVector.push_back( toPush );
        }
    }
    return true;
};

FieldOnNodesRealPtr BaseElementaryVector::assembleWithMask( const BaseDOFNumberingPtr &dofNume,
                                                            const FieldOnCellsLongPtr &maskCell,
                                                            const int &maskInve ) {

    if ( ( !dofNume ) || dofNume->isEmpty() )
        raiseAsterError( "Numerotation is empty" );

    // Create field
    FieldOnNodesRealPtr field = boost::make_shared< FieldOnNodesReal >( dofNume );

    // Elementary vector names
    std::string vectElemName = getName();
    VectorString vectElemVect( 1, vectElemName );
    char *tabNames = vectorStringAsFStrArray( vectElemVect, 19 );

    // Assembling with cell mask
    ASTERDOUBLE list_coef = 1.0;
    ASTERINTEGER typscal = 1;
    ASTERINTEGER nbElem = 1;
    std::string base( "G" );

    CALL_ASSVECWITHMASK( base.c_str(), field->getName().c_str(), &nbElem, tabNames, &list_coef,
                         dofNume->getName().c_str(), &typscal, maskCell->getName().c_str(),
                         (ASTERLOGICAL *)&maskInve );

    FreeStr( tabNames );

    return field;
};
