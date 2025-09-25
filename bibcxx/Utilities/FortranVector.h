#ifndef FORTRANVECTOR_H_
#define FORTRANVECTOR_H_

/**
 * @file FortranVector.h
 * @brief Header for Fortran tools
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

#include "astercxx.h"

/**
 * @brief add elements to set at position index of vector
 * @param pointer vector pointer
 * @param position vector index
 * @param size value size
 * @param values values to add in set
 */
void addToVectorVector( ASTERINTEGER &pointer, const ASTERINTEGER &position,
                        const ASTERINTEGER &size, const ASTERINTEGER *values );

/**
 * @brief delete set vector
 * @param pointer vector pointer to delete
 */
void deleteVectorVector( ASTERINTEGER &pointer );

/**
 * @brief create new set vector
 * @param pointer out parameter: vector pointer
 */
void getNewVectorVector( ASTERINTEGER &pointer );

/**
 * @brief get sizes of a given set vector
 * @param pointer vector pointer
 * @param vectorSize out parameter: vector size
 * @param dataSize out parameter: total element number in all sets
 */
void getVectorVectorDataSizes( ASTERINTEGER &pointer, ASTERINTEGER &vectorSize,
                               ASTERINTEGER &dataSize );

/**
 * @brief get data. WARNING: setSizes and setData must be allocated
 * @param pointer vector pointer
 * @param setSizes out parameter: vector containing all set sizes
 * @param setData out parameter: flat representation of set vector data
 */
void getVectorVectorData( ASTERINTEGER &pointer, ASTERINTEGER *setSizes, ASTERINTEGER *setData );

/**
 * @brief resize a vector
 * @param pointer vector pointer
 * @param size vector size
 */
void resizeVectorVector( ASTERINTEGER &pointer, const ASTERINTEGER &size );

#endif /* FORTRANVECTOR_H_ */
