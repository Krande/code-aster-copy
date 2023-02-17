/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org             */
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

#include "LinearAlgebra/BaseElementaryVector.h"

#include "aster_fort_calcul.h"
#include "astercxx.h"

#include "Supervis/CommandSyntax.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Tools.h"

FieldOnNodesRealPtr
BaseElementaryVector::assembleWithLoadFunctions( const BaseDOFNumberingPtr &dofNume,
                                                 const ASTERDOUBLE &time ) {
    if ( _isEmpty )
        raiseAsterError( "The ElementaryVector is empty" );

    if ( ( !dofNume ) || dofNume->isEmpty() )
        raiseAsterError( "Numerotation is empty" );

    // Elementary vector names
    std::string vectElemName = getName();

    // Pre-assembling
    std::string typres( "R" );
    std::string name( " " );
    name.resize( 24, ' ' );
    CALLO_ASASVE( vectElemName, dofNume->getName(), typres, name );

    // Get function for load
    _listOfLoads->build();
    std::string fomult( " " );
    const JeveuxVectorChar24 listOfLoadsFunc = _listOfLoads->getListOfFunctions();
    if ( listOfLoadsFunc.exists() )
        fomult = listOfLoadsFunc->getName();

    // Final assembling with load function
    FieldOnNodesRealPtr field = std::make_shared< FieldOnNodesReal >( dofNume );

    std::string detr( "D" );
    std::string base = "G";
    std::string param( "INST" );
    CALLO_ASCOVA( detr, name, fomult, param, &time, typres, field->getName(), base );

    field->updateValuePointers();
    return field;
};

FieldOnNodesRealPtr BaseElementaryVector::assembleWithMask( const BaseDOFNumberingPtr &dofNume,
                                                            const FieldOnCellsLongPtr &maskCell,
                                                            const int &maskInve ) {

    if ( ( !dofNume ) || dofNume->isEmpty() )
        raiseAsterError( "Numerotation is empty" );

    FieldOnNodesRealPtr field = std::make_shared< FieldOnNodesReal >( dofNume );

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

    field->updateValuePointers();
    return field;
};
