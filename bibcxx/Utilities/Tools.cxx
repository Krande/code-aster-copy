/**
 * @file Tools.cxx
 * @brief Implementation des outils
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

#include "Utilities/Tools.h"

#include "aster_utils.h"

std::string trim( const std::string &str, const std::string &whitespace ) {
    const std::size_t strBegin = str.find_first_not_of( whitespace );
    if ( strBegin == std::string::npos )
        return ""; // no content

    const std::size_t strEnd = str.find_last_not_of( whitespace );
    const std::size_t strRange = strEnd - strBegin + 1;

    return str.substr( strBegin, strRange );
}

std::string ljust( const std::string &str, const ASTERINTEGER &length, char fillchar ) {
    std::string tmp = str;
    tmp.resize( length, fillchar );
    return tmp;
}

std::string toUpper( const std::string &in_str ) {
    std::string str( in_str );
    std::transform( str.begin(), str.end(), str.begin(), ::toupper );
    return str;
}

std::string toLower( const std::string &in_str ) {
    std::string str( in_str );
    std::transform( str.begin(), str.end(), str.begin(), ::tolower );
    return str;
}

std::string remove_brackets( const std::string &in_str ) {
    std::string outstr;
    for ( auto ch : in_str ) {
        if ( ch == '[' ) {
            outstr += "_";
        } else if ( ch != ']' ) {
            outstr += ch;
        }
    }
    return outstr;
}

void vectorStringToFStrArray( char *tabFStr, const int flen, const VectorString &vector ) {
    int i = 0;
    for ( auto &value : vector ) {
        SetTabFStr( tabFStr, i, (char *)value.c_str(), flen );
        ++i;
    }
}

char *vectorStringAsFStrArray( const VectorString &vector, const int flen ) {
    char *tabFStr = MakeTabFStr( vector.size(), flen );
    vectorStringToFStrArray( tabFStr, flen, vector );
    return tabFStr;
}

VectorInt irange( const int begin, const int end ) {
    const int size = end - begin + 1;
    VectorInt v( size );

    int pos = 0;
    for ( int i = begin; i <= end; i++ ) {
        v[pos] = i;
        pos++;
    }

    return v;
}

VectorLong irange( const long begin, const long end ) {
    const long size = end - begin + 1;
    VectorLong v( size );

    long pos = 0;
    for ( long i = begin; i <= end; i++ ) {
        v[pos] = i;
        pos++;
    }

    return v;
}
