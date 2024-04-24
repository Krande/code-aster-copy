/**
 * @file SyntaxSaver.cxx
 * @brief Implementation de SyntaxSaver
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

#include "Utilities/SyntaxSaver.h"

#include "Utilities/SyntaxDictionary.h"
#include "pybind11/complex.h"

void SyntaxSaver::dictCopy( ListSyntaxMapContainer &toEnrich, const py::dict &toCopy ) {
    SyntaxMapContainer dict;
    for ( const auto &item : toCopy ) {
        const auto &key = item.first;
        const auto &value = item.second;
        const auto keyStr = py::cast< std::string >( key );
        if ( py::isinstance< py::str >( value ) ) {
            dict.container[keyStr] = py::cast< std::string >( value );
        } else if ( py::isinstance< py::float_ >( value ) ) {
            dict.container[keyStr] = py::cast< ASTERDOUBLE >( value );
        } else if ( py::isinstance< py::int_ >( value ) ) {
            dict.container[keyStr] = py::cast< ASTERINTEGER >( value );
        } else if ( py::isinstance< py::dict >( value ) ) {
            throw std::runtime_error( "Not allowed" );
        } else if ( py::isinstance< py::list >( value ) ) {
            auto tmp = py::cast< py::list >( value );
            auto tmp2 = tmp[0];
            if ( py::isinstance< py::str >( tmp2 ) ) {
                dict.container[keyStr] = tmp.cast< VectorString >();
            } else if ( py::isinstance< py::float_ >( tmp2 ) ) {
                dict.container[keyStr] = tmp.cast< VectorReal >();
            } else if ( py::isinstance< py::int_ >( tmp2 ) ) {
                dict.container[keyStr] = tmp.cast< VectorLong >();
            }
        } else if ( py::isinstance< py::tuple >( value ) ) {
            auto tmp0 = py::cast< py::tuple >( value );
            auto tmp = py::list( tmp0 );
            auto tmp2 = tmp[0];
            if ( py::isinstance< py::str >( tmp2 ) ) {
                dict.container[keyStr] = tmp.cast< VectorString >();
            } else if ( py::isinstance< py::float_ >( tmp2 ) ) {
                dict.container[keyStr] = tmp.cast< VectorReal >();
            } else if ( py::isinstance< py::int_ >( tmp2 ) ) {
                dict.container[keyStr] = tmp.cast< VectorLong >();
            }
        } else if ( py::isinstance< DataStructure >( value ) ) {
            auto test = py::cast< DataStructurePtr >( value );
            dict.container[keyStr] = test->getName();
        } else {
            dict.container[keyStr] = py::cast< std::complex< ASTERDOUBLE > >( value );
        }
    }
    toEnrich.push_back( dict );
}

void SyntaxSaver::listCopy( ListSyntaxMapContainer &toEnrich, const py::list &toCopy ) {
    for ( const auto &value : toCopy ) {
        auto tmp = py::cast< py::dict >( value );
        dictCopy( toEnrich, tmp );
    }
}

// SyntaxSaver::SyntaxSaver( const std::string& commandName, const ASTERINTEGER& op,
//                           py::dict syntax ): _op( op ), _syntax( commandName ) {
SyntaxSaver::SyntaxSaver( const std::string &commandName, const ASTERINTEGER &op, py::dict syntax )
    : // _op( op ), _commandName( commandName ), _keywords( syntax ) {
      _op( op ),
      _commandName( commandName ) {
    SyntaxMapContainer dict;
    for ( const auto &item : syntax ) {
        const auto &key = item.first;
        const auto &value = item.second;
        const auto keyStr = py::cast< std::string >( key );
        if ( py::isinstance< py::float_ >( value ) ) {
            dict.container[keyStr] = py::cast< ASTERDOUBLE >( value );
        } else if ( py::isinstance< py::int_ >( value ) ) {
            dict.container[keyStr] = py::cast< ASTERINTEGER >( value );
        } else if ( py::isinstance< py::str >( value ) ) {
            dict.container[keyStr] = py::cast< std::string >( value );
        } else if ( py::isinstance< py::dict >( value ) ) {
            auto tmp = py::cast< py::dict >( value );
            ListSyntaxMapContainer listeResu;
            dictCopy( listeResu, tmp );
            dict.container[keyStr] = listeResu;
        } else if ( py::isinstance< py::list >( value ) ) {
            auto tmp = py::cast< py::list >( value );
            auto tmp2 = tmp[0];
            if ( py::isinstance< py::dict >( tmp2 ) ) {
                ListSyntaxMapContainer listeResu;
                listCopy( listeResu, tmp );
                dict.container[keyStr] = listeResu;
            } else if ( py::isinstance< py::str >( tmp2 ) ) {
                dict.container[keyStr] = tmp.cast< VectorString >();
            } else if ( py::isinstance< py::float_ >( tmp2 ) ) {
                dict.container[keyStr] = tmp.cast< VectorReal >();
            } else if ( py::isinstance< py::int_ >( tmp2 ) ) {
                dict.container[keyStr] = tmp.cast< VectorLong >();
            }
        } else if ( py::isinstance< py::tuple >( value ) ) {
            auto tmp0 = py::cast< py::tuple >( value );
            auto tmp = py::list( tmp0 );
            auto tmp2 = tmp[0];
            if ( py::isinstance< py::dict >( tmp2 ) ) {
                ListSyntaxMapContainer listeResu;
                listCopy( listeResu, tmp );
                dict.container[keyStr] = listeResu;
            } else if ( py::isinstance< py::str >( tmp2 ) ) {
                dict.container[keyStr] = tmp.cast< VectorString >();
            } else if ( py::isinstance< py::float_ >( tmp2 ) ) {
                dict.container[keyStr] = tmp.cast< VectorReal >();
            } else if ( py::isinstance< py::int_ >( tmp2 ) ) {
                dict.container[keyStr] = tmp.cast< VectorLong >();
            }
        } else if ( py::isinstance< DataStructure >( value ) ) {
            auto test = py::cast< DataStructurePtr >( value );
            dict.container[keyStr] = test->getName();
        } else {
            dict.container[keyStr] = py::cast< std::complex< ASTERDOUBLE > >( value );
        }
    }
    _keywords = dict.convertToPyDict();
}
