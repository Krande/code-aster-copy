/**
 * @file Table.cxx
 * @brief Implementation de Table
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

#include "DataFields/Table.h"

/* person_in_charge: nicolas.sellenet at edf.fr */

Table::Table( const std::string &name, const std::string type )
    : DataStructure( name, 19, type ),
      _memoryLocation( JeveuxVectorChar8( getName() + ".TBBA" ) ),
      _description( JeveuxVectorLong( getName() + ".TBNP" ) ),
      _parameterDescription( JeveuxVectorChar24( getName() + ".TBLP" ) ) {};

Table::Table()
    : DataStructure( ResultNaming::getNewResultName(), 19, "TABLE" ),
      _memoryLocation( JeveuxVectorChar8( getName() + ".TBBA" ) ),
      _description( JeveuxVectorLong( getName() + ".TBNP" ) ),
      _parameterDescription( JeveuxVectorChar24( getName() + ".TBLP" ) ) {};

bool Table::build() {
    if ( _parameterDescription.exists() && _description.exists() && _parameters.size() == 0 ) {

        _parameterDescription->updateValuePointer();
        _description->updateValuePointer();
        const int nbParam = ( *_description )[0];
        _parameters.reserve( nbParam );
        for ( int i = 0; i < nbParam; ++i ) {
            std::string param = strip( ( *_parameterDescription )[i * 4].toString() );
            std::string type_name = strip( ( *_parameterDescription )[i * 4 + 1].toString() );
            if ( type_name == "R8" )
                type_name = "R";
            int typ;
            for ( typ = 0; typ < 10; typ++ ) {
                if ( JeveuxTypesNames[typ] == type_name )
                    break;
            }
            std::string name = ( *_parameterDescription )[i * 4 + 2].toString();
            std::string namev = ( *_parameterDescription )[i * 4 + 3].toString();

            if ( std::find( _parameters.begin(), _parameters.end(), param ) != _parameters.end() )
                AS_ASSERT( false )
            _parameters.push_back( param );
            _typeByString[param] = (JeveuxTypes)typ;
            switch ( typ ) {
            case JeveuxTypes::Integer:
                _columnLong[param] = JeveuxVectorLong( name );
                break;
            case JeveuxTypes::Real:
                _columnReal[param] = JeveuxVectorReal( name );
                break;
            case JeveuxTypes::Complex:
                _columnComplex[param] = JeveuxVectorComplex( name );
                break;
            case JeveuxTypes::Char8:
                _columnChar8[param] = JeveuxVectorChar8( name );
                break;
            case JeveuxTypes::Char16:
                _columnChar16[param] = JeveuxVectorChar16( name );
                break;
            case JeveuxTypes::Char24:
                _columnChar24[param] = JeveuxVectorChar24( name );
                break;
            case JeveuxTypes::Char32:
                _columnChar32[param] = JeveuxVectorChar32( name );
                break;
            case JeveuxTypes::Char80:
                _columnChar80[param] = JeveuxVectorChar80( name );
                break;
            default:
                AS_ASSERT( false );
            }
            _columnExists[param] = JeveuxVectorLong( namev );
        }
    }
#ifdef ASTER_DEBUG_CXX
    _build_called = true;
#endif
    return true;
};

int Table::getNumberOfLines() const {
    _description->updateValuePointer();
    return ( *_description )[1];
}

const VectorString Table::getParameters() const { return _parameters; };

const std::string Table::getColumnType( const std::string &param ) const {
    JeveuxTypes typ = _typeByString.at( param );
    return JeveuxTypesNames[typ];
};

const std::tuple< VectorLong, VectorLong, VectorReal, VectorComplex, VectorString >
Table::getValues( const std::string &param ) const {
    VectorLong exists;
    VectorLong vlong;
    VectorReal vreal;
    VectorComplex vcomplex;
    VectorString vstring;
#ifdef ASTER_DEBUG_CXX
    AS_ASSERT( _build_called );
#endif
    if ( std::find( _parameters.begin(), _parameters.end(), param ) != _parameters.end() ) {
        JeveuxTypes typ = _typeByString.at( param );
        _columnExists.at( param )->updateValuePointer();
        exists = _columnExists.at( param )->toVector();
        switch ( typ ) {
        case JeveuxTypes::Integer:
            _columnLong.at( param )->updateValuePointer();
            vlong = _columnLong.at( param )->toVector();
            break;
        case JeveuxTypes::Real:
            _columnReal.at( param )->updateValuePointer();
            vreal = _columnReal.at( param )->toVector();
            break;
        case JeveuxTypes::Complex:
            _columnComplex.at( param )->updateValuePointer();
            vcomplex = _columnComplex.at( param )->toVector();
            break;
        case JeveuxTypes::Char8:
            _columnChar8.at( param )->updateValuePointer();
            {
                auto v = _columnChar8.at( param )->toVector();
                vstring.reserve( v.size() );
                for ( auto c : v )
                    vstring.push_back( strip( c.toString() ) );
            }
            break;
        case JeveuxTypes::Char16:
            _columnChar16.at( param )->updateValuePointer();
            {
                auto v = _columnChar16.at( param )->toVector();
                vstring.reserve( v.size() );
                for ( auto c : v )
                    vstring.push_back( strip( c.toString() ) );
            }
            break;
        case JeveuxTypes::Char24:
            _columnChar24.at( param )->updateValuePointer();
            {
                auto v = _columnChar24.at( param )->toVector();
                vstring.reserve( v.size() );
                for ( auto c : v )
                    vstring.push_back( strip( c.toString() ) );
            }
            break;
        case JeveuxTypes::Char32:
            _columnChar32.at( param )->updateValuePointer();
            {
                auto v = _columnChar32.at( param )->toVector();
                vstring.reserve( v.size() );
                for ( auto c : v )
                    vstring.push_back( strip( c.toString() ) );
            }
            break;
        case JeveuxTypes::Char80:
            _columnChar80.at( param )->updateValuePointer();
            {
                auto v = _columnChar80.at( param )->toVector();
                vstring.reserve( v.size() );
                for ( auto c : v )
                    vstring.push_back( strip( c.toString() ) );
            }
            break;
        default:
            AS_ASSERT( false );
        }
    }
    return std::tuple( exists, vlong, vreal, vcomplex, vstring );
}

Table::~Table() {
#ifdef ASTER_DEBUG_CXX
    std::cout << "DEBUG: Table.destr: " << this->getName() << std::endl;
#endif
};

TableOfFunctions::TableOfFunctions( const std::string &name ) : Table( name, "TABLE_FONCTION" ) {};

TableOfFunctions::TableOfFunctions()
    : Table( ResultNaming::getNewResultName(), "TABLE_FONCTION" ) {};

void TableOfFunctions::addFunction( GenericFunctionPtr func ) {
    _vecOfFunctions.push_back( func );
};

GenericFunctionPtr TableOfFunctions::getFunction( int pos ) const {
    if ( pos < int( _vecOfFunctions.size() ) )
        return _vecOfFunctions[pos];
    return BaseFunctionPtr( nullptr );
};

int TableOfFunctions::getNumberOfFunctions() const { return _vecOfFunctions.size(); };
