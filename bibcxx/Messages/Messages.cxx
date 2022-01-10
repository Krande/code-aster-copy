/**
 * @file Message.cxx
 * @brief Fichier entete de la class Messages
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

#include "Messages/Messages.h"
#include "aster_fort_utils.h"

void UTMESS( char *error, char *message ) { CALL_UTMESS( error, message ); }

void UTMESS( const char *error, const char *message ) { UTMESS( (char *)error, (char *)message ); }

void UTMESS( const std::string &error, const std::string &message ) {
    UTMESS( &error[0], &message[0] );
}
