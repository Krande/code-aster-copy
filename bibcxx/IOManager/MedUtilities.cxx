/**
 * @file MedProfile.h
 * @brief Fichier entete de la classe MedProfile
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "MedUtilities.h"

#include "Utilities/Tools.h"

#include <cstring>
#include <iostream>

// aslint: disable=C3012

std::vector< std::string > splitChar( char *toSplit, int nbElem, int size ) {
    std::vector< std::string > toReturn;
    for ( int i = 0; i < nbElem; ++i ) {
        char *tmp = toSplit + i * size;
        const auto str = std::string( tmp, std::min( strlen( tmp ), (size_t)size ) );
        toReturn.push_back( strip( str ) );
    }
    return toReturn;
};

std::pair< int, int > splitEntitySet( int nbElemT, int rank, int nbProcs ) {
    int nbElemL = nbElemT / nbProcs;
    int start = rank * nbElemL + 1;
    const auto end = nbProcs * nbElemL;
#ifdef ASTER_DEBUG_CXX
    int lastNbElemL = nbElemL - ( end - nbElemT );
    int verif = nbElemL * ( nbProcs - 1 ) + lastNbElemL;
    if ( verif != nbElemT ) {
        const auto str1 = std::to_string( verif );
        const auto str2 = std::to_string( nbElemT );
        throw std::runtime_error( "Error in splitting " + str1 + " " + str2 );
    }
#endif
    if ( rank == nbProcs - 1 ) {
        nbElemL = nbElemL - ( end - nbElemT );
    }
    return { nbElemL, start };
};
