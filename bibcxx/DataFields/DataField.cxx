/**
 * @file DataField.cxx
 * @brief Implementation de DataField vide car DataField est un template
 * @author Nicolas Pignet
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

#include "DataFields/DataField.h"

#include "aster_fort_utils.h"

std::string DataField::getFieldType() const {
    const std::string questi1( "TYPE_CHAMP" );
    const std::string typeco( "CHAMP" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "F" );

    CALLO_DISMOI( questi1, getName(), typeco, &repi, repk, arret, &ier );

    return strip( repk.toString() );
};

std::string DataField::getFieldScalar() const {
    const std::string questi1( "TYPE_SCA" );
    const std::string typeco( "CHAMP" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "F" );

    CALLO_DISMOI( questi1, getName(), typeco, &repi, repk, arret, &ier );

    return strip( repk.toString() );
};
