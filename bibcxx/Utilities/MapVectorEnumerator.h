#ifndef MAPVECTORENUMERATOR_H_
#define MAPVECTORENUMERATOR_H_

/**
 * @file MapVectorEnumerator.h
 * @brief Iterate on vector like on map (with index, value)
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

#include <iostream>
#include <map>
#include <utility>
#include <vector>

/**
 * @brief enumerate function to iterate on vector with (index, value) non const case
 * @author Nicolas Sellenet
 */
template < typename T >
auto enumerate( std::vector< T > &vec ) {
    struct iterator {
        size_t index;
        typename std::vector< T >::iterator it;
        std::pair< size_t, T & > pair; // Stockage de la paire

        iterator( size_t i, typename std::vector< T >::iterator iter )
            : index( i ), it( iter ), pair( i, *iter ) {}

        bool operator!=( const iterator &other ) const { return it != other.it; }
        void operator++() {
            ++index;
            ++it;
            if ( it != typename std::vector< T >::iterator {} ) {
                pair = { index, *it };
            }
        }
        void operator--() {
            --index;
            --it;
            if ( it != typename std::vector< T >::iterator {} ) {
                pair = { index, *it };
            }
        }
        auto &operator*() { return pair; }
    };

    struct wrapper {
        std::vector< T > &vec;
        auto begin() { return iterator { 0, vec.begin() }; }
        auto end() { return iterator { vec.size(), vec.end() }; }
    };

    return wrapper { vec };
}

/**
 * @brief enumerate function to iterate on vector with (index, value) const case
 * @author Nicolas Sellenet
 */
template < typename T >
auto enumerate( const std::vector< T > &vec ) {
    struct iterator {
        size_t index;
        typename std::vector< T >::const_iterator it;
        std::pair< size_t, const T & > pair;

        iterator( size_t i, typename std::vector< T >::const_iterator iter )
            : index( i ), it( iter ), pair( i, *iter ) {}

        bool operator!=( const iterator &other ) const { return it != other.it; }
        void operator++() {
            ++index;
            ++it;
            if ( it != typename std::vector< T >::const_iterator {} ) {
                pair.first = index;
                // NOTE, object still const, only pointer change !
                const_cast< T & >( pair.second ) = *it;
            }
        }
        void operator--() {
            --index;
            --it;
            if ( it != typename std::vector< T >::const_iterator {} ) {
                pair.first = index;
                // NOTE, object still const, only pointer change !
                const_cast< T & >( pair.second ) = *it;
            }
        }
        auto &operator*() const { return pair; }
    };

    struct wrapper {
        const std::vector< T > &vec;
        auto begin() const { return iterator { 0, vec.begin() }; }
        auto end() const { return iterator { vec.size(), vec.end() }; }
    };

    return wrapper { vec };
}

/**
 * @brief enumerate function to iterate on map with (key, value) non const case
 * @author Nicolas Sellenet
 */
template < typename K, typename V >
auto enumerate( std::map< K, V > &m ) {
    struct iterator {
        typename std::map< K, V >::iterator it;
        bool operator!=( const iterator &other ) const { return it != other.it; }
        void operator++() { ++it; }
        void operator--() { --it; }
        auto &operator*() { return *it; }
    };

    struct wrapper {
        std::map< K, V > &m;
        auto begin() { return iterator { m.begin() }; }
        auto end() { return iterator { m.end() }; }
    };

    return wrapper { m };
}

/**
 * @brief enumerate function to iterate on map with (key, value) const case
 * @author Nicolas Sellenet
 */
template < typename K, typename V >
auto enumerate( const std::map< K, V > &m ) {
    struct iterator {
        typename std::map< K, V >::const_iterator it;
        bool operator!=( const iterator &other ) const { return it != other.it; }
        void operator++() { ++it; }
        void operator--() { --it; }
        auto &operator*() const { return *it; }
    };

    struct wrapper {
        const std::map< K, V > &m;
        auto begin() const { return iterator { m.begin() }; }
        auto end() const { return iterator { m.end() }; }
    };

    return wrapper { m };
}

#endif /* MAPVECTORENUMERATOR_H_ */
