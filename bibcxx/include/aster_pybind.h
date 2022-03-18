#pragma once

/**
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

#include "astercxx.h"

#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"

#include <pybind11/complex.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

/**
 * @brief Factory for '__init__' constructor without 'DSTypePtr'.
 */
template < typename DSType, typename... Args >
static std::shared_ptr< DSType > initFactoryPtr( Args... args ) {
    return std::shared_ptr< DSType >( std::make_shared< DSType >( args... ) );
};

namespace pybind11 {
namespace detail {

/**
 * @brief Converter for JeveuxVector.
 */
template < typename T >
struct type_caster< JeveuxVector< T > > {
  public:
    PYBIND11_TYPE_CASTER( JeveuxVector< T >, const_name( "JeveuxVector" ) );

    bool load( handle /* src */, bool ) { return false; }

    static handle cast( const JeveuxVector< T > &vect, return_value_policy /* policy */,
                        handle /* parent */ ) {
        py::list pylist;
        vect->updateValuePointer();
        for ( int i = 0; i < vect->size(); ++i ) {
            pylist.append( ( *vect )[i] );
        }
        return pylist.inc_ref();
    }
};

/**
 * @brief Converter for JeveuxCollection.
 */
template < typename T >
struct type_caster< JeveuxCollection< T > > {
  public:
    PYBIND11_TYPE_CASTER( JeveuxCollection< T >, const_name( "JeveuxCollection" ) );

    bool load( handle /* src */, bool ) { return false; }

    static handle cast( const JeveuxCollection< T > &coll, return_value_policy /* policy */,
                        handle /* parent */ ) {
        py::list pylist;
        if ( !coll->build() || coll->size() < 0 ) {
            return pylist.inc_ref();
        }
        for ( int i = 0; i < coll->size(); ++i ) {
            const auto &obj = coll->getObject( i + 1 );
            py::list items;
            for ( int j = 0; j < obj.size(); ++j ) {
                items.append( obj[j] );
            }
            items.inc_ref();
            pylist.append( items );
        }
        return pylist.inc_ref();
    }
};

} // namespace detail
} // namespace pybind11
