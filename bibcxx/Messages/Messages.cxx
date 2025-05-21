/**
 * @file Message.cxx
 * @brief Fichier entete de la class Messages
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

#include "Messages/Messages.h"

#include "aster_fort_utils.h"

void UTMESS( const std::string &typm, const std::string &idmess ) {
    if ( typm == "A" || typm == "I" ) {
        CALL_UTMESS( (char *)typm.c_str(), (char *)idmess.c_str() );
    } else {
        raiseAsterError( idmess );
    }
}
void UTMESS( char *typm, char *idmess ) { UTMESS( std::string( typm ), std::string( idmess ) ); }

void UTMESS( const char *typm, const char *idmess ) { UTMESS( (char *)typm, (char *)idmess ); }

void UtmessCore( const std::string &typm, const std::string &idmess, const VectorString &vec ) {
    ASTERINTEGER n0 = 0, n1 = 1, ibid = 0;
    ASTERDOUBLE rbid = 0.;
    std::string typm2( typm ), idmess2( idmess );
    char *valk, *fname;
    fname = MakeBlankFStr( 1 );
    valk = MakeTabFStr( vec.size(), VALK_SIZE );
    for ( int i = 0; i < vec.size(); ++i ) {
        SetTabFStr( valk, i, vec[i].data(), VALK_SIZE );
    }
    CALL_UTMESS_CORE( typm2.data(), idmess2.data(), &n1, valk, &n0, &ibid, &n0, &rbid, &n0, fname );
    FreeStr( valk );
    FreeStr( fname );
}
