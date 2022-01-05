/**
 * @file ContactNew.cxx
 * @brief Implementation de Contact
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

#include "Contact/ContactNew.h"
#include "Messages/Messages.h"

ContactNew::ContactNew( const std::string name, const ModelPtr model )
    : DataStructure( name, 8, "CHAR_CONT" ), _model( model ),
      _FEDesc( boost::make_shared< FiniteElementDescriptor >( getName() + ".CHME.LIGRE",
                                                              _model->getMesh() ) ),
      _verbosity( 1 ), _friction( false ), _smoothing( false ) {
    // model has to be mechanics
    if ( !_model->isMechanical() )
        UTMESS( "F", "CONTACT1_2" );
};

bool ContactNew::build() {
    // build FiniteElementDescriptor
    return true;    
}
