/**
 * @file NonLinearResult.cxx
 * @brief Implementation de NonLinearResult
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

#include "Results/NonLinearResult.h"
#include "Supervis/Exceptions.h"

void NonLinearResult::setContact( const ContactPtr contact,
                                         const ASTERINTEGER& rank ) {
    if( !contact )
      raiseAsterError( "ValueError: Contact is empty" );

    _mapContact[rank] = contact;
    const auto fed = contact->getFiniteElementDescriptor();
    _fieldBuidler.addFiniteElementDescriptor( fed );
};

void NonLinearResult::setContact( const ContactPtr contact ) {
    _serialNumber->updateValuePointer();
    ASTERINTEGER nbRanks = _serialNumber->size();
    for ( ASTERINTEGER rank = 0; rank < nbRanks; ++rank ) {
        const ASTERINTEGER iordr = ( *_serialNumber )[rank];
        if ( _mapContact.find( iordr ) == _mapContact.end() )
            setContact( contact, iordr );
    }
};
