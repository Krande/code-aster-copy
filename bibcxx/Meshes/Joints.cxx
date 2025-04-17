/**
 * @file ParallelDOFNumbering.cxx
 * @brief Implementation de ParallelDOFNumbering
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "Meshes/Joints.h"

#ifdef ASTER_HAVE_MPI

#include "aster_mpi.h"

#include "ParallelUtilities/AsterMPI.h"
#include "Supervis/Exceptions.h"

#include <algorithm>

Joints::Joints() : Joints( DataStructureNaming::getNewName() ) {};

Joints::Joints( const std::string name )
    : DataStructure( name, 19, "DOMJOINTS" ),
      _layer( JeveuxVectorLong( getName() + ".NBLG" ) ),
      _domj( JeveuxVectorLong( getName() + ".DOMJ" ) ),
      _send( JeveuxCollectionLong( getName() + ".SEND" ) ),
      _recv( JeveuxCollectionLong( getName() + ".RECV" ) ),
      _procGroupIds( JeveuxVectorShort( getName() + ".PGID" ) ),
      _groupComm( JeveuxVectorLong( getName() + ".GCOM" ) ) {
    buildGroup();
};

void Joints::buildGroup() {
    _procGroupIds->deallocate();
    _groupComm->deallocate();
    const auto aster_comm = aster_get_comm_world();
    const MPI_Comm comm = aster_comm->id;
    int rank, size;
    aster_get_mpi_info( aster_comm, &rank, &size );
    VectorInt toGather, vecOut, idsInParent, groupProcIds;
    if ( _domj->exists() )
        toGather.push_back( rank );
    else
        toGather.push_back( -1 );
    AsterMPI::all_gather( toGather, vecOut );
    int cmpt = 0;
    for ( auto id : vecOut ) {
        if ( id != -1 ) {
            idsInParent.push_back( cmpt );
            groupProcIds.push_back( id );
            ++cmpt;
        } else {
            idsInParent.push_back( -1 );
        }
    }
    *_procGroupIds = idsInParent;

    if ( _procGroupIds->exists() ) {
        _jointsMPIGroup = MPIGroupPtr( new MPIGroup() );
        _jointsMPIGroup->buildFromProcsVector( comm, groupProcIds );
        _groupComm->allocate( 1 );
        const auto tmp = MPI_Comm_c2f( _jointsMPIGroup->getCommunicator() );
        ( *_groupComm )[0] = tmp;
    }
};

void Joints::setOppositeDomains( const VectorLong &oppositeDomains ) {
    if ( oppositeDomains.size() != 0 )
        ( *_domj ) = oppositeDomains;
    buildGroup();
};

const JeveuxVectorLong &Joints::getOppositeDomains() const {
    if ( _domj->exists() ) {
        _domj->updateValuePointer();
    }

    return _domj;
}

bool Joints::build() {
    _send->build();
    _recv->build();

    return true;
}

void Joints::setSendedElements( const VectorOfVectorsLong &send ) {
    if ( send.size() > 0 ) {
        _send->allocate( send.size() );

        ASTERINTEGER i = 1;
        for ( auto &send_i : send ) {
            if ( send_i.size() > 0 ) {
                auto obj = _send->allocateObject( i, send_i );
            }
            i++;
        }
    }
};

void Joints::setReceivedElements( const VectorOfVectorsLong &recv ) {
    if ( recv.size() > 0 ) {
        _recv->allocate( recv.size() );

        ASTERINTEGER i = 1;
        for ( auto &recv_i : recv ) {
            if ( recv_i.size() > 0 ) {
                auto obj = _recv->allocateObject( i, recv_i );
            }
            i++;
        }
    }
};

VectorLong Joints::getSendedElements( const ASTERINTEGER &id ) const {
    const auto &obj = ( *_send )[id + 1];
    obj->updateValuePointer();
    return obj->toVector();
};

VectorLong Joints::getReceivedElements( const ASTERINTEGER &id ) const {
    const auto &obj = ( *_recv )[id + 1];
    obj->updateValuePointer();
    return obj->toVector();
};

void Joints::setNumberOfGhostLayer( const ASTERINTEGER &nb_layer ) {
    if ( !_layer->exists() ) {
        _layer->allocate( 1 );
    }
    _layer->updateValuePointer();

    if ( nb_layer < 1 ) {
        raiseAsterError( "ValueError: ghost layer number must be at least 1." );
    }

    ( *_layer )[0] = nb_layer;
};

ASTERINTEGER Joints::getNumberOfGhostLayer( void ) const {
    _layer->updateValuePointer();
    return ( *_layer )[0];
};

#endif /* ASTER_HAVE_MPI */
