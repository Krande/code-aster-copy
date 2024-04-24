/**
 * @file ObjectBalancer.cxx
 * @brief Implementation of an object balancer
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

#include "ParallelUtilities/ObjectBalancer.h"

#ifdef ASTER_HAVE_MPI

void ObjectBalancer::prepareCommunications() {
    if ( !_sendDefined ) {
        throw std::runtime_error( "The definition of elementary sends must finished"
                                  " end before calling prepareCommunications" );
    }
    _graph->synchronizeOverProcesses();
    const auto rank = getMPIRank();
    // Communicate what to send and what to receive
    for ( const auto [tag, proc] : *_graph ) {
        if ( proc == -1 )
            continue;
        VectorInt tmp( 1, -1 );
        if ( rank > proc ) {
            tmp[0] = _sendList[proc].size();
            AsterMPI::send( tmp, proc, tag );
            AsterMPI::receive( tmp, proc, tag );
#ifdef ASTER_DEBUG_CXX
            std::cout << "#" << rank << " received " << tmp[0] << " from #" << proc << std::endl;
#endif
            _recvSize[proc] = tmp[0];
        } else {
            AsterMPI::receive( tmp, proc, tag );
#ifdef ASTER_DEBUG_CXX
            std::cout << "#" << rank << " received " << tmp[0] << " from #" << proc << std::endl;
#endif
            _recvSize[proc] = tmp[0];
            tmp[0] = _sendList[proc].size();
            AsterMPI::send( tmp, proc, tag );
        }
    }
    std::set< int > delInterKeep;
    std::set_intersection( _toDelete.begin(), _toDelete.end(), _toKeep.begin(), _toKeep.end(),
                           std::inserter( delInterKeep, delInterKeep.begin() ) );
    if ( delInterKeep.size() != 0 ) {
        throw std::runtime_error( "Inconsistent communication definition."
                                  " Some elements are to keep and to delete at the same time." );
    }
    std::set< int > sendInterKeep;
    std::set_intersection( _toSend.begin(), _toSend.end(), _toKeep.begin(), _toKeep.end(),
                           std::inserter( sendInterKeep, sendInterKeep.begin() ) );
    std::set< int > sendDiffDel;
    std::set_difference( _toSend.begin(), _toSend.end(), _toDelete.begin(), _toDelete.end(),
                         std::inserter( sendDiffDel, sendDiffDel.begin() ) );
    const auto nbProcs = getMPISize();

    // Compute size delta for vectors
    _sizeDelta = 0;
    for ( int iProc = 0; iProc < nbProcs; ++iProc ) {
        _sizeDelta += _recvSize[iProc];
    }
    _sizeDelta -= ( sendDiffDel.size() + _toDelete.size() - sendInterKeep.size() );
    _isOk = true;
};

void ObjectBalancer::balanceObjectOverProcesses( const MeshCoordinatesFieldPtr &coordsIn,
                                                 MeshCoordinatesFieldPtr &coordsOut ) const {
    if ( !_isOk )
        throw std::runtime_error( "ObjectBalancer not prepared" );
    auto valuesIn = coordsIn->getValues();
    const auto vecSize = valuesIn->size();

    auto valuesOut = coordsOut->getValues();
    valuesOut->allocate( vecSize + 3 * _sizeDelta );
    valuesOut->updateValuePointer();
    balanceSimpleVectorOverProcesses< ASTERDOUBLE, 3 >( valuesIn->getDataPtr(), vecSize,
                                                        valuesOut->getDataPtr() );
    coordsOut->buildDescriptor();
};

#endif /* ASTER_HAVE_MPI */
