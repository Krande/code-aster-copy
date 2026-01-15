#ifndef MEDUTILITIES_H_
#define MEDUTILITIES_H_

/**
 * @file MedUtilities.h
 * @brief Fichier entete de la classe MedUtilities
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
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

#include "astercxx.h"

#ifdef ASTER_HAVE_MED
#include "med.h"

#include <string>
#include <vector>

static const VectorInt asterTypeList = { 1,  2,  4,  6,  7,  9,  11, 12, 14, 16,
                                         18, 19, 20, 21, 22, 23, 24, 25, 26, 27 };

static const std::map< int, med_geometry_type > asterMedMatching = {
    { 1, 1 },    { 2, 102 },  { 4, 103 },  { 6, 104 },  { 7, 203 },  { 9, 206 },  { 11, 207 },
    { 12, 204 }, { 14, 208 }, { 16, 209 }, { 18, 304 }, { 19, 310 }, { 20, 306 }, { 21, 315 },
    { 22, 318 }, { 23, 305 }, { 24, 313 }, { 25, 308 }, { 26, 320 }, { 27, 327 }
};

const std::set< med_int > medTypeToRenumber = { 304, 308, 305, 306, 310, 320, 313, 315, 318, 327 };

// aslint: disable=C3012

/** @brief split char* in nbElem std::string of size size */
std::vector< std::string > splitChar( char *toSplit, int nbElem, int size );

/** @brief copy VectorString to fixed size char* */
char *stringVectorToChar( const VectorString &vec, int fixedSize );

/** @brief parallel split a set of element (used for med filters) */
std::pair< int, int > splitEntitySet( int nbElemT, int rank, int nbProcs );

template < std::size_t N, const int indices[N] >
void applyPermutation( const med_int *in, med_int *out ) {
    for ( size_t i = 0; i < N; i++ ) {
        out[i] = in[indices[i]];
    }
};
#endif

#endif /* MEDUTILITIES_H_ */
