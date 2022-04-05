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

bool ContactZone::build() {
    auto mesh = getMesh();

    // check that there is no nodes in common
    auto slaveNodes_lc = mesh->getNodesFromCells( getSlaveGroupOfCells(), false, true );
    auto masterNodes_lc = mesh->getNodesFromCells( getMasterGroupOfCells(), false, true );

    VectorLong commonNodes;

    if ( mesh->isParallel() ) {
#ifdef ASTER_HAVE_MPI
        VectorLong slaveNodes_gl;
        AsterMPI::all_gather( slaveNodes_lc, slaveNodes_gl );
        commonNodes = set_intersection( slaveNodes_gl, masterNodes_lc );
#endif
    } else {
        commonNodes = set_intersection( slaveNodes_lc, masterNodes_lc );
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

    // Update  master nodes
    _masterNodes = std::move( masterNodes_lc );

    // check mesh orientation (normals)
    if ( checkNormals() ) {
        std::string slave = ljust( getSlaveGroupOfCells(), 24, ' ' );
        std::string master = ljust( getMasterGroupOfCells(), 24, ' ' );
        CALL_CHECKNORMALS( _model->getName().c_str(), slave.c_str(), master.c_str() );
    }

    // Be sure that the lists of cells are updated
    updateMasterCells();
    updateSlaveCells();

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
    masterCells.reserve(_masterCells.size());
    for (auto cell: _masterCells)
        masterCells.push_back( cell + 1 );
    CALL_CNCINV( getMesh()->getName().c_str(), masterCells.data(), &nbMaster, base.c_str(),
                 _masterInverseConnectivity->getName().c_str() );
    _masterInverseConnectivity->build();

    // create slave inverse connectivity
    ASTERINTEGER nbSlave = getSlaveCells().size();

    VectorLong slaveCells;
    slaveCells.reserve(_slaveCells.size());
    for (auto cell: _slaveCells)
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
        masterCells.reserve(_masterCells.size());
        for (auto cell: _masterCells)
            masterCells.push_back( cell + 1 );
        CALL_CNVOIS( getMesh()->getName(), masterCells.data(), invmcn_name, &nbMaster, &ind_min,
                     &ind_max, mn_name );

        _masterNeighbors->build();
    }

    // get slave neighbors

    ASTERINTEGER nbSlave = getSlaveCells().size();
    if ( nbSlave > 0 ) {
        ind_max = *std::max_element( _slaveCells.begin(), _slaveCells.end() ) + 1;
        ind_min = *std::min_element( _slaveCells.begin(), _slaveCells.end() ) + 1;

        std::string invscn_name = ljust( _slaveInverseConnectivity->getName(), 24, ' ' );
        std::string sn_name = ljust( _slaveNeighbors->getName(), 24, ' ' );

        VectorLong slaveCells;
        slaveCells.reserve(_slaveCells.size());
        for (auto cell: _slaveCells)
            slaveCells.push_back( cell + 1 );
        CALL_CNVOIS( getMesh()->getName(), slaveCells.data(), invscn_name, &nbSlave, &ind_min,
                     &ind_max, sn_name );

        _slaveNeighbors->build();
    }

    return true;
}
