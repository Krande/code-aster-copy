/**
 * @file ContactZone.cxx
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

#include "Contact/ContactZone.h"

#include "aster_fort_mesh.h"

#include "Messages/Messages.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

ContactZone::ContactZone( const std::string name, const ModelPtr model )
    : DataStructure( name, 8, "CHAR_CONT_ZONE" ),
      _model( model ),
      _verbosity( 1 ),
      _checkNormal( true ),
      _smoothing( false ),
      _contParam( std::make_shared< ContactParameter >() ),
      _fricParam( std::make_shared< FrictionParameter >() ),
      _pairParam( std::make_shared< PairingParameter >() ),
      _masterInverseConnectivity( JeveuxCollectionLong( getName() + ".CM" ) ),
      _slaveInverseConnectivity( JeveuxCollectionLong( getName() + ".CS" ) ),
      _masterNeighbors( JeveuxCollectionLong( getName() + ".MN" ) ),
      _slaveNeighbors( JeveuxCollectionLong( getName() + ".SN" ) ) {
    // model has to be mechanics
    if ( !_model->isMechanical() )
        UTMESS( "F", "CONTACT1_2" );
};

void ContactZone::setSlaveGroupOfCells( const std::string &slave ) {
    if ( getMesh()->hasGroupOfCells( slave ) ) {
        _slaveNodes = getMesh()->getNodesFromCells( slave, false, true );
        _slaveCells = getMesh()->getCells( slave );
        _slaveGrp = slave;
    } else {
        throw std::runtime_error( "The given group " + slave + " doesn't exist in mesh" );
    }
};

void ContactZone::setMasterGroupOfCells( const std::string &master ) {
    if ( getMesh()->hasGroupOfCells( master ) ) {
        _masterNodes = getMesh()->getNodesFromCells( master, false, true );
        _masterCells = getMesh()->getCells( master );
        _masterGrp = master;
    } else {
        throw std::runtime_error( "The given group " + master + " doesn't exist in mesh" );
    }
};

void ContactZone::setExcludedSlaveGroupOfCells( const VectorString &excluded_slave ) {
    auto mesh = getMesh();
    for ( auto &name : excluded_slave ) {
        if ( !( getMesh()->hasGroupOfCells( name ) ) ) {
            throw std::runtime_error( "The group " + name + " doesn't exist in mesh" );
        }

        auto it = _slaveCellsExcluded.end();
        VectorLong sans_gr_i = mesh->getCells( name );
        _slaveCellsExcluded.insert( it, sans_gr_i.begin(), sans_gr_i.end() );
    }
    _slaveCells = set_intersection( _slaveCells, _slaveCellsExcluded );
    _slaveNodes = getMesh()->getNodesFromCells( _slaveCells );
};

bool ContactZone::build() {
    auto mesh = getMesh();

    // check that there is no nodes in common
    VectorLong commonNodes;

    if ( mesh->isParallel() ) {
#ifdef ASTER_HAVE_MPI
        VectorLong slaveNodes_gl;
        AsterMPI::all_gather( _slaveNodes, slaveNodes_gl );
        commonNodes = set_intersection( slaveNodes_gl, _masterNodes );
#endif
    } else {
        commonNodes = set_intersection( _slaveNodes, _masterNodes );
    }

    ASTERINTEGER size_inter_gl = commonNodes.size();

    // share error
#ifdef ASTER_HAVE_MPI
    if ( mesh->isParallel() ) {
        ASTERINTEGER size_inter_lc = size_inter_gl;
        size_inter_gl = 0;
        AsterMPI::all_reduce( size_inter_lc, size_inter_gl, MPI_SUM );
    }
#endif

    if ( size_inter_gl > 0 ) {
        UTMESS( "F", "CONTACT1_1" );
    }

    // check mesh orientation (normals)
    if ( checkNormals() ) {
        CALL_CHECKNORMALS( _model->getName().c_str(), ljust( _slaveGrp, 24, ' ' ).c_str(),
                           ljust( _masterGrp, 24, ' ' ).c_str() );
    }

    // build inverse connvectivity
    buildInverseConnectivity();

    // build master and slave  Cells Neighbors
    buildCellsNeighbors();

    return true;
}

ASTERBOOL ContactZone::buildInverseConnectivity() {
    // create master inverse connectivity
    ASTERINTEGER nbMaster = getMasterCells().size();
    std::string base( "G" );

    VectorLong masterCells;
    masterCells.reserve( _masterCells.size() );
    for ( auto cell : _masterCells )
        masterCells.push_back( cell + 1 );
    CALL_CNCINV( getMesh()->getName().c_str(), masterCells.data(), &nbMaster, base.c_str(),
                 _masterInverseConnectivity->getName().c_str() );
    _masterInverseConnectivity->build();

    // create slave inverse connectivity
    ASTERINTEGER nbSlave = _slaveCells.size();

    VectorLong slaveCells;
    slaveCells.reserve( _slaveCells.size() );
    for ( auto cell : _slaveCells )
        slaveCells.push_back( cell + 1 );
    CALL_CNCINV( getMesh()->getName().c_str(), slaveCells.data(), &nbSlave, base.c_str(),
                 _slaveInverseConnectivity->getName().c_str() );
    _slaveInverseConnectivity->build();
    return true;
}

ASTERBOOL ContactZone::buildCellsNeighbors() {

    ASTERINTEGER ind_max, ind_min;
    // get master neighbors
    ASTERINTEGER nbMaster = getMasterCells().size();
    if ( nbMaster > 0 ) {
        ind_max = *std::max_element( _masterCells.begin(), _masterCells.end() ) + 1;
        ind_min = *std::min_element( _masterCells.begin(), _masterCells.end() ) + 1;

        std::string invmcn_name = ljust( _masterInverseConnectivity->getName(), 24, ' ' );
        std::string mn_name = ljust( _masterNeighbors->getName(), 24, ' ' );

        VectorLong masterCells;
        masterCells.reserve( _masterCells.size() );
        for ( auto cell : _masterCells )
            masterCells.push_back( cell + 1 );
        CALL_CNVOIS( getMesh()->getName(), masterCells.data(), invmcn_name, &nbMaster, &ind_min,
                     &ind_max, mn_name );

        _masterNeighbors->build();
    }

    // get slave neighbors

    ASTERINTEGER nbSlave = _slaveCells.size();
    if ( nbSlave > 0 ) {
        ind_max = *std::max_element( _slaveCells.begin(), _slaveCells.end() ) + 1;
        ind_min = *std::min_element( _slaveCells.begin(), _slaveCells.end() ) + 1;

        std::string invscn_name = ljust( _slaveInverseConnectivity->getName(), 24, ' ' );
        std::string sn_name = ljust( _slaveNeighbors->getName(), 24, ' ' );

        VectorLong slaveCells;
        slaveCells.reserve( _slaveCells.size() );
        for ( auto cell : _slaveCells )
            slaveCells.push_back( cell + 1 );
        CALL_CNVOIS( getMesh()->getName(), slaveCells.data(), invscn_name, &nbSlave, &ind_min,
                     &ind_max, sn_name );

        _slaveNeighbors->build();
    }

    return true;
}

VectorLong ContactZone::getMasterCellsFromNode( const ASTERINTEGER &i ) const {
    auto vct = ( *_masterInverseConnectivity )[i + 1].toVector();
    std::transform( vct.begin(), vct.end(), vct.begin(), [this]( ASTERINTEGER k ) -> ASTERINTEGER {
        return k > 0 ? _masterCells[k - 1] : 0;
    } );
    return vct;
}

VectorLong ContactZone::getSlaveCellsFromNode( const ASTERINTEGER &i ) const {
    auto vct = ( *_slaveInverseConnectivity )[i + 1].toVector();
    std::transform( vct.begin(), vct.end(), vct.begin(), [this]( ASTERINTEGER k ) -> ASTERINTEGER {
        return k > 0 ? _slaveCells[k - 1] : 0;
    } );
    return vct;
}

VectorLong ContactZone::getMasterCellNeighbors( const ASTERINTEGER &i ) const {
    ASTERINTEGER ind_min = *std::min_element( _masterCells.begin(), _masterCells.end() );
    ASTERINTEGER ind_max = *std::max_element( _masterCells.begin(), _masterCells.end() );

    if ( i < ind_min || i > ind_max )
        throw std::out_of_range( " the master cell's number should be"
                                 " between " +
                                 std::to_string( ind_min ) + " and " + std::to_string( ind_max ) );

    auto vct = ( *_masterNeighbors )[i - ind_min + 1].toVector();
    vct.erase( std::remove_if( vct.begin(), vct.end(), []( ASTERINTEGER &i ) { return i == 0; } ),
               vct.end() );
    std::transform( vct.begin(), vct.end(), vct.begin(),
                    []( ASTERINTEGER k ) -> ASTERINTEGER { return k - 1; } );
    return vct;
}

VectorLong ContactZone::getSlaveCellNeighbors( const ASTERINTEGER &i ) const {

    ASTERINTEGER ind_min = *std::min_element( _slaveCells.begin(), _slaveCells.end() );
    ASTERINTEGER ind_max = *std::max_element( _slaveCells.begin(), _slaveCells.end() );

    if ( i < ind_min || i > ind_max )
        throw std::out_of_range( " the slave cell's number should be"
                                 " between " +
                                 std::to_string( ind_min ) + " and " + std::to_string( ind_max ) );

    auto vct = ( *_slaveNeighbors )[i - ind_min + 1].toVector();
    vct.erase( std::remove_if( vct.begin(), vct.end(), []( ASTERINTEGER &i ) { return i == 0; } ),
               vct.end() );
    std::transform( vct.begin(), vct.end(), vct.begin(),
                    []( ASTERINTEGER k ) -> ASTERINTEGER { return k - 1; } );
    return vct;
}
