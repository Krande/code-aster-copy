/**
 * @file DataStructure.cxx
 * @brief Implementation des fonctions membres de DataStructure
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

#include "DataStructures/DataStructure.h"

#include "shared_vars.h"

#include "MemoryManager/JeveuxString.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

#include <algorithm>

// Constructor
DataStructure::DataStructure( const std::string name, const int nameLength, const std::string type )
    : _name( name ), _sdj( py::none() ), _cache( py::none() ) {
    _name.resize( nameLength, ' ' );

    // #ifdef ASTER_DEBUG_CXX_OBJECTS
    //     std::cout << "Creating " << strip( this->getName() )
    //             << " <" << type << "> "
    //             << this->getUserName() << std::endl;
    // #endif

    std::string name19( _name );
    name19.resize( 19, ' ' );
    _tco = JeveuxVectorChar24( name19 + "._TCO" );
    if ( !_tco.exists() && name != "" ) {
        _tco->allocate( 1 );
        if ( type.size() <= 8 && type != "FORMULE" )
            ( *_tco )[0] = std::string( strip( type ) + "_SDASTER" );
        else
            ( *_tco )[0] = type;
    }

    // title
    _title = JeveuxVectorChar80( name19 + ".TITR" );
}

// Constructor
DataStructure::DataStructure( const int nameLength, const std::string type )
    : DataStructure::DataStructure( DataStructureNaming::getNewName( nameLength ), nameLength,
                                    type ) {}

// Copy Constructor
DataStructure::DataStructure( const DataStructure &other )
    : DataStructure::DataStructure( DataStructureNaming::getNewName( other.getName().size() ),
                                    other.getName().size(), other.getType() ) {
    _user_name = other._user_name;
    *( _title ) = *( other._title );
    _depsVector = other._depsVector;
}

// Move Constructor
DataStructure::DataStructure( DataStructure &&other ) {
    auto name = other.getName();
    other._name = "";
    _name = name;
    _user_name = other._user_name;
    _tco = other._tco;
    _title = other._title;
    _depsVector = other._depsVector;
    _sdj = other._sdj;
}

// Assignment operator
DataStructure &DataStructure::operator=( DataStructure &&other ) {
    if ( this != &other ) {
        auto name = other.getName();
        other._name = "";
        _name = name;
        _user_name = other._user_name;
        _tco = other._tco;
        _title = other._title;
        _depsVector = other._depsVector;
        _sdj = other._sdj;
    }
    return *this;
}

// Destructor
DataStructure::~DataStructure() {
    // user/main datastructures aka 'concept' (== without ".").
    bool mainDs = _name.find( "." ) == std::string::npos && _name[0] != '&';
    std::string nameWithoutBlanks = strip( _name );
    // empty name or no memory manager : skip silently
    if ( nameWithoutBlanks == "" || get_sh_jeveux_status() != 1 ) {
        return;
    }
    // Allow to see when the datastructure is really deleted.
    // In case of embraced datastructures, '_tco' is deallocated the first time
    // (no type)
#ifdef ASTER_DEBUG_CXX_OBJECTS
    if ( mainDs && this->getType() != "not_found" ) {
        // Too low-level to call UTMESS.
        std::cout << "Deleting " << strip( this->getName() ) << " <" << this->getType() << "> "
                  << this->getUserName() << std::endl;
    }
#endif
    // Destruction
    _tco->deallocate();
    _title->deallocate();
#ifdef ASTER_DEBUG_CXX_OBJECTS
    std::string base( " " );
    ASTERINTEGER pos = 1;
    ASTERINTEGER nbval2 = 0;
    ASTERINTEGER retour = 0;
    JeveuxChar24 nothing( " " );
    if ( nameWithoutBlanks == "&2" ) {
        retour = 1;
    }
    CALLO_JELSTC( base, nameWithoutBlanks, &pos, &nbval2, nothing, &retour );
    if ( retour != 0 ) {
        JeveuxVectorChar24 test( "&&TMP" );
        test->allocate( -retour );
        ASTERINTEGER nbval2 = -retour;
        CALLO_JELSTC( base, nameWithoutBlanks, &pos, &nbval2, ( *test )[0], &retour );
        std::cout << "DEBUG: Remaining jeveux objects in " << _name << std::endl;
        std::cout << "DEBUG: List of objects:" << std::endl;
        for ( int i = 0; i < retour; ++i )
            std::cout << "DEBUG:   - " << ( *test )[i].toString() << std::endl;
    }
#endif
};

void DataStructure::addDependency( const DataStructurePtr &ds ) {
    int idx;
    int size( _depsVector.size() );
    for ( idx = 0; idx < size; idx++ ) {
        if ( ds == _depsVector[idx] ) {
            break;
        }
    }
    if ( idx == size ) {
        _depsVector.push_back( ds );
    }
}

void DataStructure::removeDependency( const DataStructurePtr &ds ) {
    int idx;
    int size( _depsVector.size() );
    for ( idx = 0; idx < size; idx++ ) {
        if ( ds == _depsVector[idx] ) {
            break;
        }
    }
    if ( idx < size ) {
        _depsVector.erase( _depsVector.begin() + idx );
    }
}

void DataStructure::resetDependencies() { _depsVector.clear(); }

std::vector< DataStructure::DataStructurePtr > DataStructure::getDependencies() const {
    return _depsVector;
}

void DataStructure::debugPrint( int logicalUnit, bool synchro ) const {
    ASTERINTEGER unit, niveau, ipos, True, False;
    unit = (ASTERINTEGER)logicalUnit;
    niveau = 2;
    True = 1;
    False = 0;
    ipos = 1;
    JeveuxString< 1 > base( " " );
    JeveuxString< 3 > no( "NON" );
    std::string nameWithoutBlanks = strip( _name );
#ifdef ASTER_HAVE_MPI
    int rank = getMPIRank();
    int nbProcs = getMPISize();
    if ( !synchro ) {
        rank = 0;
        nbProcs = 1;
    }
    for ( int iProc = 0; iProc < nbProcs; ++iProc ) {
        if ( iProc == rank ) {
#endif /* ASTER_HAVE_MPI */
            try {
                CALLO_UTIMSD( &unit, &niveau, &False, &True, nameWithoutBlanks, &ipos, base, no );
            } catch ( ... ) {
                throw std::runtime_error( "debugPrint failed!" );
            }
#ifdef ASTER_HAVE_MPI
        }
        if ( synchro )
            AsterMPI::barrier();
    }
#endif /* ASTER_HAVE_MPI */
};

void DataStructure::setType( const std::string newType ) {
    CALL_JEMARQ();
    _tco->updateValuePointer();
    if ( newType.size() <= 8 && newType != "FORMULE" )
        ( *_tco )[0] = std::string( strip( newType ) + "_SDASTER" );
    else
        ( *_tco )[0] = newType;
    CALL_JEDEMA();
};

void DataStructure::setUserName( const std::string name ) { _user_name = name; }

void DataStructure::setTitle( const std::string title ) {
    CALL_JEMARQ();
    if ( !_title.exists() )
        _title->allocate( 1 );

    _title->updateValuePointer();

    ( *_title )[0] = title;
    CALL_JEDEMA();
}

std::string DataStructure::getTitle() {
    if ( !_title.exists() ) {
        return std::string();
    }

    CALL_JEMARQ();
    _title->updateValuePointer();

    std::string title = ( *_title )[0].toString();
    CALL_JEDEMA();

    return title;
}
