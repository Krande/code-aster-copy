/**
 * @file ContactPairing.cxx
 * @brief Implementation de Contact
 * @section LICENCE
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

#include "Contact/ContactPairing.h"

ContactPairing::ContactPairing( const std::string name, const std::vector< ContactZonePtr > zones,
                                const BaseMeshPtr mesh )
    : DataStructure( name, 8, "PAIRING_SD" ), _zones( zones ), _mesh( mesh ) {
    if ( !_mesh || _mesh->isParallel() )
        raiseAsterError( "Mesh is empty or is parallel " );

    _newCoordinates = std::make_shared< MeshCoordinatesField >( *( _mesh->getCoordinates() ) );

    // be sure that zones is not empty and get size of zones and resize
    if ( _zones.empty() )
        raiseAsterError( "ContactZone vector is empty " );
    int size_zones = zones.size();

    // resize pairing quantities
    _nbPairs.resize( size_zones );
    _listOfPairs.resize( size_zones );
    _nbIntersectionPoints.resize( size_zones );
    _slaveIntersectionPoints.resize( size_zones );
};

ASTERBOOL ContactPairing::computeZone( ASTERINTEGER i ) {

    if ( i < 0 || i >= _zones.size() ) {
        throw std::out_of_range( "The zone index should be between 0  and " +
                                 std::to_string( _zones.size() - 1 ) );
    }

    // get and define some input parameters
    VectorLong eleMaster = _zones[i]->getMasterCells();
    VectorLong NodesMaster = _zones[i]->getMasterNodes();
    VectorLong eleSlave = _zones[i]->getSlaveCells();
    ASTERINTEGER nbCellMaster = eleMaster.size();
    ASTERINTEGER nbNodeMaster = NodesMaster.size();
    ASTERINTEGER nbCellSlave = eleSlave.size();
    std::string pair_method;

    // update the numbering for fortran
    std::for_each( eleMaster.begin(), eleMaster.end(), []( ASTERINTEGER &d ) { d += 1; } );
    std::for_each( NodesMaster.begin(), NodesMaster.end(), []( ASTERINTEGER &d ) { d += 1; } );
    std::for_each( eleSlave.begin(), eleSlave.end(), []( ASTERINTEGER &d ) { d += 1; } );

    // get pairing method
    auto variant = _zones[i]->getPairingParameter()->getAlgorithm();
    if ( variant == PairingAlgo::Mortar ) {
        pair_method = ljust( "ROBUSTE", 24, ' ' );
    } else {
        AS_ABORT( "Not expected" );
    }

    // tolerence
    ASTERDOUBLE pair_tole = 1e-8;

    // set pairs numbers to 0
    ASTERINTEGER nb_pairs = 0;

    // output paramaters as C pointers
    ASTERINTEGER *pairs = NULL;
    ASTERINTEGER *nbInterPoints = NULL;
    ASTERDOUBLE *InterSlavePoints = NULL;

    CALLO_APLCPGN( _mesh->getName(), _newCoordinates->getName(), _zones[i]->getName(), pair_method,
                   &pair_tole, &nbCellMaster, eleMaster.data(), &nbCellSlave, eleSlave.data(),
                   NodesMaster.data(), &nbNodeMaster, &nb_pairs, &pairs, &nbInterPoints,
                   &InterSlavePoints );

    // clearZone
    this->clearZone( i );

    // fill the pairing quantities
    _nbPairs[i] = nb_pairs;
    _listOfPairs[i] = VectorLong( pairs, pairs + 2 * nb_pairs );
    _nbIntersectionPoints[i] = VectorLong( nbInterPoints, nbInterPoints + nb_pairs );
    _slaveIntersectionPoints[i] = VectorReal( InterSlavePoints, InterSlavePoints + 16 * nb_pairs );

    // update numerotation

    std::transform( _listOfPairs[i].begin(), _listOfPairs[i].end(), _listOfPairs[i].begin(),
                    []( ASTERINTEGER &i ) -> ASTERINTEGER { return --i; } );

    // free temporary quantities
    if ( nb_pairs > 0 ) {
        free( pairs );
        free( nbInterPoints );
        free( InterSlavePoints );
    }

    return true;
}

ASTERBOOL ContactPairing::compute() {
    for ( int i = 0; i < _zones.size(); i++ ) {
        computeZone( i );
    }

    return true;
}

void ContactPairing::clearZone( ASTERINTEGER i ) {

    // swap is recommended to release memory
    _listOfPairs[i].clear();
    _nbIntersectionPoints[i].clear();
    _slaveIntersectionPoints[i].clear();

    _nbPairs.at( i ) = 0;
}

std::vector< std::pair< ASTERINTEGER, ASTERINTEGER > >
ContactPairing::getListOfPairsOfZone( ASTERINTEGER zone_index ) const {

    std::vector< std::pair< ASTERINTEGER, ASTERINTEGER > > tmp;
    ASTERINTEGER nbPairs = getNumberOfPairsOfZone( zone_index );
    tmp.reserve( nbPairs );

    for ( auto i = 0; i < nbPairs; i++ ) {
        tmp.push_back( std::make_pair( _listOfPairs[zone_index][2 * i],
                                       _listOfPairs[zone_index][2 * i + 1] ) );
    }
    return tmp;
}

std::vector< VectorReal >
ContactPairing::getSlaveIntersectionPoints( ASTERINTEGER zone_index ) const {

    std::vector< VectorReal > ret;
    ASTERINTEGER nbPairs = getNumberOfPairsOfZone( zone_index );
    ret.reserve( nbPairs );

    auto iter = _slaveIntersectionPoints[zone_index].begin();
    for ( auto i = 0; i < nbPairs; i++ ) {
        ret.push_back( VectorReal( iter + 16 * i, iter + 16 * ( i + 1 ) ) );
    }

    return ret;
}
