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
    _masterIntersectionPoints.resize( size_zones );
    _quadraturePoints.resize( size_zones );
};

ASTERBOOL ContactPairing::compute( ASTERINTEGER i ) {

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

    // update the numbering
    std::for_each( eleMaster.begin(), eleMaster.end(), []( ASTERINTEGER &d ) { d += 1; } );
    std::for_each( NodesMaster.begin(), NodesMaster.end(), []( ASTERINTEGER &d ) { d += 1; } );
    std::for_each( eleSlave.begin(), eleSlave.end(), []( ASTERINTEGER &d ) { d += 1; } );

    // get pairing method
    ContactVariant variant = _zones[i]->getContactParameter()->getVariant();
    if ( variant == ContactVariant::Robust ) {
        pair_method = ljust( "ROBUSTE", 24, ' ' );
    } else if ( variant == ContactVariant::Rapide ) {
        pair_method = ljust( "RAPIDE", 24, ' ' );
    } else {
        pair_method = ljust( "ROBUSTE", 24, ' ' );
    }

    // tolerence
    ASTERDOUBLE pair_tole = 0.001;

    // set pairs numbers to 0
    ASTERINTEGER nb_pairs = 0;

    // output paramaters as C pointers
    ASTERINTEGER *pairs = NULL;
    ASTERINTEGER *nbInterPoints = NULL;
    ASTERDOUBLE *InterSlavePoints = NULL;
    ASTERDOUBLE *InterMasterPoints = NULL;
    ASTERDOUBLE *gaussPoints = NULL;

    CALLO_APLCPGN( _mesh->getName(), _newCoordinates->getName(), _zones[i]->getName(), pair_method,
                   &pair_tole, &nbCellMaster, eleMaster.data(), &nbCellSlave, eleSlave.data(),
                   NodesMaster.data(), &nbNodeMaster, &nb_pairs, &pairs, &nbInterPoints,
                   &InterSlavePoints, &InterMasterPoints, &gaussPoints );

    // clearZone
    this->clearZone( i );

    // fill the pairing quantities
    _nbPairs[i] = nb_pairs;
    _listOfPairs[i] = VectorLong( pairs, pairs + 2 * nb_pairs );
    _nbIntersectionPoints[i] = VectorLong( nbInterPoints, nbInterPoints + nb_pairs );
    _slaveIntersectionPoints[i] = VectorReal( InterSlavePoints, InterSlavePoints + 16 * nb_pairs );
    _masterIntersectionPoints[i] =
        VectorReal( InterMasterPoints, InterMasterPoints + 16 * nb_pairs );
    _quadraturePoints[i] = VectorReal( gaussPoints, gaussPoints + 72 * nb_pairs );

    // update numerotation

    std::transform( _listOfPairs[i].begin(), _listOfPairs[i].end(), _listOfPairs[i].begin(),
                    []( ASTERINTEGER &i ) -> ASTERINTEGER { return --i; } );

    // free temporary quantities
    free( pairs );
    free( nbInterPoints );
    free( InterSlavePoints );
    free( InterMasterPoints );
    free( gaussPoints );

    return true;
}

ASTERBOOL ContactPairing::clearZone( ASTERINTEGER i ) {

    // swap is recommended to release memory
    VectorLong().swap( _listOfPairs[i] );
    VectorLong().swap( _nbIntersectionPoints[i] );
    VectorReal().swap( _slaveIntersectionPoints[i] );
    VectorReal().swap( _masterIntersectionPoints[i] );
    VectorReal().swap( _quadraturePoints[i] );

    _nbPairs.at( i ) = 0;

    return true;
}
