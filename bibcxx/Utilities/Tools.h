#ifndef TOOLS_H_
#define TOOLS_H_

/**
 * @file Tools.h
 * @brief Fichier entete des outils
 * @author Nicolas Sellenet
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

#include <string>
#include <set>

#include "astercxx.h"

std::string trim( const std::string &str, const std::string &whitespace = " \t" );

/**
 * @brief irange Create a vector of integer from begin to end (included).
 *      for exemple {-1, 0, 1, 2, 3}
 */
VectorInt irange(const int begin, const int end);

VectorLong irange(const long begin, const long end);


/**
 * @brief vectorStringToFStr Create an array of Fortran strings from a vector of strings.
 *      The output array must be freed by the caller.
 */
char *vectorStringAsFStrArray( const VectorString &vector, const int size );


// Set and sort a vector
template<typename T>
std::vector<T>
unique(const std::vector<T>& vec)
{
    // make unique & sort
    std::set<T> s;
    std::copy(vec.begin(), vec.end(), std::inserter(s, s.end()));

    // recopy
    std::vector<T> r;
    r.resize(s.size());
    std::copy(s.begin(), s.end(), r.begin());

    return r;
}

// Get unique list of Aster Concept in a map (indexed by rank)
template<typename T>
std::vector<T>
unique(const std::map< int, T >& _map)
{
    std::map< std::string, T> unique_map;
    for(auto it : _map)
    {
        unique_map[it.second->getName()] = it.second;
    }

    std::vector< T > ret;
    ret.reserve(unique_map.size());
    for(auto it : unique_map)
    {
        ret.push_back(it.second);
    }

    return ret;
}

#endif /* TOOLS_H_ */
