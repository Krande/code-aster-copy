/**
 * @file FortranTools.cxx
 * @brief Implementation des outils
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

#include "Utilities/FortranVector.h"

typedef std::vector< VectorInt > VectorOfVectorsInt;

/** @brief vector of all set vectors */
std::vector< VectorOfVectorsInt > vectorVector;

/** @brief map from pointer to vector index in vectorVector */
std::map< ASTERINTEGER, ASTERINTEGER > ref;

void addToVectorVector( ASTERINTEGER &pointer, const ASTERINTEGER &position,
                        const ASTERINTEGER &size, const ASTERINTEGER4 *values ) {
    auto ptr = (VectorOfVectorsInt *)pointer;
    auto &vec = *ptr;
    auto &curSet = vec[position];
    for ( int i = 0; i < size; ++i ) {
        curSet.push_back( values[i] );
    }
}

void deleteVectorVector( ASTERINTEGER &pointer ) {
    auto pos = ref[pointer];
    vectorVector[pos] = VectorOfVectorsInt();
    ref.erase( ref.find( pointer ) );
}

void getNewVectorVector( ASTERINTEGER &pointer ) {
    vectorVector.push_back( VectorOfVectorsInt() );
    auto *tmp = &vectorVector[vectorVector.size() - 1];
    pointer = (ASTERINTEGER)tmp;
    ref[pointer] = vectorVector.size() - 1;
}

void getVectorVectorData( ASTERINTEGER &pointer, ASTERINTEGER *setSizes, ASTERINTEGER4 *setData ) {
    auto ptr = (VectorOfVectorsInt *)pointer;
    auto &vec = *ptr;
    int cmpt1 = 0, cmpt2 = 0, cumSize = 1;
    for ( const auto &curSet : vec ) {
        setSizes[cmpt1] = cumSize;
        cumSize += curSet.size();
        for ( const auto &val : curSet ) {
            setData[cmpt2] = val;
            ++cmpt2;
        }
        ++cmpt1;
    }
}

void getVectorVectorDataSizes( ASTERINTEGER &pointer, ASTERINTEGER &vectorSize,
                               ASTERINTEGER &dataSize ) {
    auto ptr = (VectorOfVectorsInt *)pointer;
    auto &vec = *ptr;
    vectorSize = vec.size();
    dataSize = 0;
    for ( const auto &curSet : vec ) {
        dataSize += curSet.size();
    }
}

void resizeVectorVector( ASTERINTEGER &pointer, const ASTERINTEGER &size ) {
    auto ptr = (VectorOfVectorsInt *)pointer;
    for ( int i = 0; i < size; ++i ) {
        ptr->push_back( VectorInt() );
    }
}

void sortUniqueVectorVector( ASTERINTEGER &pointer ) {
    auto ptr = (VectorOfVectorsInt *)pointer;
    auto &vec = *ptr;
    int cmpt1 = 0, cmpt2 = 0, cumSize = 1;
    for ( auto &curSet : vec ) {
        sort( curSet.begin(), curSet.end() );

        // Group unique elements together
        auto it = unique( curSet.begin(), curSet.end() );

        // Erase duplicates
        curSet.erase( it, curSet.end() );
    }
}

extern "C" void DEFPPPP( VECTOR_VECTOR_ADD_VALUES, vector_vector_add_values, ASTERINTEGER *ptr,
                         const ASTERINTEGER *position, const ASTERINTEGER *size,
                         const ASTERINTEGER4 *values ) {
    addToVectorVector( *ptr, *position, *size, values );
}

extern "C" void DEFP( VECTOR_VECTOR_DELETE, vector_vector_delete, ASTERINTEGER *ptr ) {
    deleteVectorVector( *ptr );
}

extern "C" void DEFP( VECTOR_VECTOR_NEW, vector_vector_new, ASTERINTEGER *ptr ) {
    getNewVectorVector( *ptr );
}

extern "C" void DEFPPP( VECTOR_VECTOR_GET_DATA, vector_vector_get_data, ASTERINTEGER *ptr,
                        ASTERINTEGER *setSizes, ASTERINTEGER4 *setData ) {
    getVectorVectorData( *ptr, setSizes, setData );
}

extern "C" void DEFPPP( VECTOR_VECTOR_GET_DATA_SIZES, vector_vector_get_data_sizes,
                        ASTERINTEGER *ptr, ASTERINTEGER *vectorSize, ASTERINTEGER *dataSize ) {
    getVectorVectorDataSizes( *ptr, *vectorSize, *dataSize );
}

extern "C" void DEFPP( VECTOR_VECTOR_RESIZE, vector_vector_resize, ASTERINTEGER *ptr,
                       const ASTERINTEGER *size ) {
    resizeVectorVector( *ptr, *size );
}

extern "C" void DEFP( VECTOR_VECTOR_SORT_UNIQUE, vector_vector_sort_unique, ASTERINTEGER *ptr ) {
    sortUniqueVectorVector( *ptr );
}
