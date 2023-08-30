/**
 * @file ObjectBalancer.cxx
 * @brief Implementation of an object balancer
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
    const auto tmp = std::set< int >( _toSend.begin(), _toSend.end() );
    std::set< int > intersect;
    std::set_intersection( tmp.begin(), tmp.end(), _toKeep.begin(), _toKeep.end(),
                           std::inserter( intersect, intersect.begin() ) );
    const auto nbProcs = getMPISize();

    // Compute size delta for vectors
    _sizeDelta = 0;
    for ( int iProc = 0; iProc < nbProcs; ++iProc ) {
        _sizeDelta += _recvSize[iProc];
    }
    _sizeDelta -= ( _toSend.size() - intersect.size() );
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

#ifdef ASTER_HAVE_MED
MedVectorPtr
ObjectBalancer::balanceMedVectorOverProcessesWithRenumbering( const MedVectorPtr &vecIn ) const {
    MedVectorPtr vecOut( new MedVector() );
    balanceObjectOverProcesses3( *vecIn, *vecOut, DummyMaskDouble() );
    if ( _renumbering.size() == 0 ) {
        return vecOut;
    }
    const auto size = vecOut->size();
    MedVectorPtr vecOut2( new MedVector() );
    vecOut2->setComponentNumber( vecOut->getComponentNumber() );
    vecOut2->setSize( size );
    if ( _renumbering.size() != size )
        throw std::runtime_error( "Sizes not matching" );
    for ( int i = 0; i < size; ++i ) {
        const auto newId = _renumbering[i] - 1;
        vecOut2->setElement( newId, vecOut->getElement( i ) );
    }
    vecOut2->endDefinition();
    for ( int i = 0; i < size; ++i ) {
        const auto newId = _renumbering[i] - 1;
        const auto nbCmp = vecOut->getElement( i );
        const auto &elInR = ( *vecOut )[i];
        auto &elOutR = ( *vecOut2 )[newId];
        for ( int j = 0; j < nbCmp; ++j ) {
            elOutR[j] = elInR[j];
        }
    }
    return vecOut2;
};
#endif

#endif /* ASTER_HAVE_MPI */
