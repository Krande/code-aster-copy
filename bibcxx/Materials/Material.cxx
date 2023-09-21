/**
 * @file Material.cxx
 * @brief Implementation of material datastructure.
 * @section LICENCE
 * Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
 * This file is part of code_aster.
 *
 * code_aster is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * code_aster is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "astercxx.h"

#include "Materials/Material.h"

#include "aster_fort_material.h"

#include "Utilities/Tools.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

/* Lists */
MaterialListReal::MaterialListReal( const std::string &name, const MaterialListReal &toCopy )
    : MaterialListReal( name, toCopy._list->toVector() ) {}

MaterialListFunc::MaterialListFunc( const std::string &name, const VectorString &vect )
    : DataStructure( name, 16, "MATER_PROP_LISTF" ),
      _list( JeveuxVectorChar8( getName() + ".LISV_FO", vect.size() ) ) {
    int idx = 0;
    for ( auto elt : vect ) {
        ( *_list )[idx] = elt;
        ++idx;
    }
}

MaterialListFunc::MaterialListFunc( const std::string &name, const MaterialListFunc &toCopy )
    : DataStructure( name, 16, "MATER_PROP_LISTF" ),
      _list( JeveuxVectorChar8( getName() + ".LISV_FO", toCopy._list->size() ) ) {
    int idx = 0;
    for ( auto elt : toCopy._list->toVector() ) {
        ( *_list )[idx] = elt;
        ++idx;
    }
}

/* Material Property */
MaterialProperties::MaterialProperties( const std::string &name )
    : DataStructure( name, 19, "MATER_PROP" ),
      _valR( JeveuxVectorReal( getName() + ".VALR" ) ),
      _valC( JeveuxVectorComplex( getName() + ".VALC" ) ),
      _valK( JeveuxVectorChar16( getName() + ".VALK" ) ) {}

MaterialProperties::MaterialProperties( const std::string &name, const MaterialProperties &toCopy )
    : MaterialProperties( name ) {
    *( _valR ) = *( toCopy._valR );
    *( _valC ) = *( toCopy._valC );
    *( _valK ) = *( toCopy._valK );
    if ( toCopy._ordr.exists() ) {
        _ordr = JeveuxVectorChar16( getName() + ".ORDR" );
        *( _ordr ) = *( toCopy._ordr );
    }
    if ( toCopy._kord.exists() ) {
        _kord = JeveuxVectorLong( getName() + ".KORD" );
        *( _kord ) = *( toCopy._kord );
    }
    if ( getValueString( 0 ) == "LISTE_COEF" ) {
        std::string listName = getValueString( 1 );
        listName.replace( 0, 8, toCopy.getName() );
    }
}

MaterialProperties::MaterialProperties( const std::string &name, const int nbParam,
                                        const VectorReal valR, const VectorComplex valC,
                                        const VectorString valK, const VectorString ordr,
                                        const VectorLong kord )
    : DataStructure( name, 19, "MATER_PROP" ),
      _valR( JeveuxVectorReal( getName() + ".VALR", nbParam ) ),
      _valC( JeveuxVectorComplex( getName() + ".VALC", nbParam ) ),
      _valK( JeveuxVectorChar16( getName() + ".VALK", 2 * nbParam ) ) {
    ( *_valR ) = 0.;
    _valR->setSize( valR.size() );
    for ( int idx = 0; idx < valR.size(); idx++ ) {
        ( *_valR )[idx] = valR[idx];
    }

    ( *_valC ) = std::complex( 0., 0. );
    _valC->setSize( valC.size() );
    for ( int idx = 0; idx < valC.size(); idx++ ) {
        ( *_valC )[valR.size() + idx] = valC[idx];
    }

    _valK->setSize( valK.size() );
    int idx = 0;
    for ( auto elt : valK ) {
        ( *_valK )[idx] = elt;
        ++idx;
    }

    if ( ordr.size() > 0 ) {
        _ordr = JeveuxVectorChar16( getName() + ".ORDR", ordr.size() );
        int idx = 0;
        for ( auto elt : ordr ) {
            ( *_ordr )[idx] = elt;
            ++idx;
        }
        _kord = JeveuxVectorLong( getName() + ".KORD", kord );
    }
}

ASTERDOUBLE MaterialProperties::getValueReal( int idx ) {
    _valR->updateValuePointer();
    return ( *_valR )[idx];
}

ASTERCOMPLEX MaterialProperties::getValueComplex( int idx ) {
    _valC->updateValuePointer();
    return ( *_valC )[idx];
}

std::string MaterialProperties::getValueString( int idx ) {
    _valK->updateValuePointer();
    return ( *_valK )[idx];
}

/* Material */
Material::Material( const std::string &name )
    : DataStructure( name, 8, "MATER" ),
      _names( JeveuxVectorChar32( getName() + ".MATERIAU.NOMRC" ) ),
      _rdep( std::make_shared< Function >( getName() + ".&&RDEP" ) ) {
    build();
}

Material::Material( const Material &toCopy ) : Material( ResultNaming::getNewResultName() ) {
    *( _names ) = *( toCopy._names );
    int rdepSize = toCopy._rdep->maximumSize();
    if ( rdepSize > 0 ) {
        _rdep->allocate( rdepSize );
        _rdep->setParameterName( "EPSI" );
        _rdep->setResultName( toCopy._rdep->getResultName() );
    }

    int nbMat = size();
    for ( int i = 0; i < nbMat; i++ ) {
        auto prop = *( toCopy._prop[i] );
        auto copy = std::make_shared< MaterialProperties >( _cptName( i + 1 ), prop );
        _prop.push_back( copy );
    }

    for ( auto &ds : toCopy.getDependencies() ) {
        addDependency( ds );
    }
}

Material::Material( const Material &toCopy, const VectorString propIgnored )
    : Material( ResultNaming::getNewResultName() ) {

    int rdepSize = toCopy._rdep->maximumSize();
    if ( rdepSize > 0 ) {
        _rdep->allocate( rdepSize );
        _rdep->setParameterName( "EPSI" );
        _rdep->setResultName( toCopy._rdep->getResultName() );
    }

    auto &tocopynames = *( toCopy._names );
    int nbMat = tocopynames.size();
    int nbCpt = 1;
    for ( int i = 0; i < nbMat; i++ ) {

        if ( find( propIgnored.begin(), propIgnored.end(), tocopynames[i].rstrip() ) ==
             propIgnored.end() ) {
            _names->push_back( tocopynames[i] );
            auto prop = *( toCopy._prop[i] );
            auto copy = std::make_shared< MaterialProperties >( _cptName( nbCpt ), prop );
            _prop.push_back( copy );
            nbCpt++;
        }
    }

    for ( auto &ds : toCopy.getDependencies() ) {
        addDependency( ds );
    }
}

bool Material::build() {
    int nbMat = size();
    if ( nbMat > 0 ) {
        auto first = _names->front().rstrip();
        if ( first == "ELAS_COQMU" || first == "THER_COQMU" ) {
            nbMat = 1;
        }
    }
    for ( int i = 0; i < nbMat; i++ ) {
        auto prop = std::make_shared< MaterialProperties >( _cptName( i + 1 ) );
        _prop.push_back( prop );
        if ( prop->getValueString( 0 ) == "LISTE_COEF" ) {
            _nameList.push_back( prop->getValueString( 1 ) );
        }
    }
    return true;
}

ASTERINTEGER Material::size() {
    if ( !_names.exists() ) {
        return 0;
    }
    _names->updateValuePointer();
    return _names->size();
}

std::string Material::_storeListReal( VectorReal vect ) {
    std::ostringstream numb;
    numb << std::setw( 7 ) << std::setfill( '0' ) << _lisvR.size() + _lisvF.size() + 1;
    std::string name = getName() + "." + numb.str();

    _lisvR.push_back( std::make_shared< MaterialListReal >( name, vect ) );
    _nameList.push_back( name );
    return name;
}

std::string Material::_storeListFunc( VectorString vect ) {
    std::ostringstream numb;
    numb << std::setw( 7 ) << std::setfill( '0' ) << _lisvR.size() + _lisvF.size() + 1;
    std::string name = getName() + "." + numb.str();

    _lisvF.push_back( std::make_shared< MaterialListFunc >( name, vect ) );
    _nameList.push_back( name );
    return name;
}

void Material::_addProperties( const std::string name, int nbParam, VectorReal valR,
                               VectorComplex valC, VectorString valK, VectorString ordr,
                               VectorLong kord ) {
    if ( _names.exists() && std::find( _names.begin(), _names.end(), name ) != _names.end() ) {
        throw std::runtime_error( "Properties for '" + name + "' are already defined" );
    }
    _names->push_back( name );
    std::string cpt = _cptName( _names->size() );
    _prop.push_back(
        std::make_shared< MaterialProperties >( cpt, nbParam, valR, valC, valK, ordr, kord ) );
    if ( name == "ELAS_ISTR" || name == "ELAS_ORTH" ) {
        CALLO_CHECK_ANISO( name, cpt );
    }
}

void Material::_setTractionFunction( std::string name, std::string keyword,
                                     GenericFunctionPtr &trac ) {
    ASTERINTEGER dummy = 0;
    CALLO_RCSTOC_VERIF( trac->getName(), keyword, name, &dummy );
    _rdep->allocate( trac->maximumSize() );
    _rdep->setParameterName( "EPSI" );
    _rdep->setResultName( trac->getResultName() );
}

std::string Material::getListName( const int position ) {
    std::string name = "";
    if ( position >= int( _prop.size() ) )
        // throw std::runtime_error( "getListName: Out of bound" );
        return name;
    auto prop = _prop[position];
    if ( prop->getValueString( 0 ) == "LISTE_COEF" &&
         std::find( _nameList.begin(), _nameList.end(), prop->getValueString( 1 ) ) !=
             _nameList.end() ) {
        name = prop->getValueString( 1 );
    }
    return name;
}

VectorString Material::getMaterialNames() {
    VectorString names;
    names.reserve( size() );
    for ( auto &elt : *_names ) {
        names.push_back( strip( elt.toString() ) );
    }
    return names;
}

ASTERDOUBLE Material::getValueReal( const std::string matName, const std::string propName ) {
    // get object name
    auto prop = matByName( matName );
    if ( !prop ) {
        throw std::runtime_error( "material not found: " + matName );
    }
    std::string objName = "";
    int nbObj = prop->getNumberOfReal();
    for ( int i = 0; i < nbObj; i++ ) {
        if ( strip( prop->getValueString( i ) ) == propName ) {
            return prop->getValueReal( i );
        }
    }
    throw std::runtime_error( "property not found: " + propName );
};

ASTERCOMPLEX Material::getValueComplex( const std::string matName, const std::string propName ) {
    // get object name
    auto prop = matByName( matName );
    if ( !prop ) {
        throw std::runtime_error( "material not found: " + matName );
    }
    std::string objName = "";
    int first = prop->getNumberOfReal();
    int nbObj = prop->getNumberOfComplex();
    for ( int i = 0; i < nbObj; i++ ) {
        if ( strip( prop->getValueString( first + i ) ) == propName ) {
            return prop->getValueComplex( first + i );
        }
    }
    throw std::runtime_error( "property not found: " + propName );
};

GenericFunctionPtr Material::getFunction( const std::string matName, const std::string propName ) {
    // get object name
    auto prop = matByName( matName );
    if ( !prop ) {
        return nullptr;
    }
    std::string objName = "";
    int first = prop->getNumberOfReal() + prop->getNumberOfComplex();
    int nbObj = prop->getNumberOfObjects();
#ifdef ASTER_DEBUG_CXX
    std::cout << "DEBUG: search from pos: " << first << " with nbObj: " << nbObj << std::endl;
#endif
    for ( int i = 0; i < nbObj; i++ ) {
        if ( strip( prop->getValueString( first + i ) ) == propName ) {
            objName = prop->getValueString( first + nbObj + i );
            break;
        }
    }
    if ( objName == "" ) {
#ifdef ASTER_DEBUG_CXX
        std::cout << "DEBUG: property not found: " + propName << std::endl;
#endif
        return nullptr;
    }
    // getDependencyByName
    GenericFunctionPtr ds( nullptr );
    for ( auto &elt : getDependencies() ) {
        if ( strip( elt->getName() ) == strip( objName ) ) {
            ds = std::static_pointer_cast< GenericFunction >( elt );
            break;
        }
    }
    if ( !ds ) {
#ifdef ASTER_DEBUG_CXX
        std::cout << "DEBUG: object '" + objName + "' not found in dependencies" << std::endl;
#endif
        return nullptr;
    }
    return ds;
};

/* Material, private functions */
std::string Material::_cptName( int idx ) {
    std::ostringstream numb;
    numb << std::setw( 6 ) << std::setfill( '0' ) << idx;
    return getName() + ".CPT." + numb.str();
}

MaterialPropertiesPtr Material::matByName( std::string matName ) {
    // private: do not use it with *_COQMU
    int idx = -1, i = 0;
    for ( auto &elt : *_names ) {
        if ( elt.rstrip() == matName ) {
            idx = i;
            break;
        }
        i++;
    }
    if ( idx < 0 ) {
        return nullptr;
    }
    return _prop[idx];
}
