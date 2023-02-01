#ifndef OBJECTBALANCER_H_
#define OBJECTBALANCER_H_

/**
 * @file ObjectBalancer.h
 * @brief Header of an object balancer
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
#include "aster_mpi.h"
#include "astercxx.h"

#include "MemoryManager/JeveuxString.h"
#include "MemoryManager/JeveuxVector.h"
#include "ParallelUtilities/CommGraph.h"

/**
 * @class CommGraph
 * @brief Class used to "balance" objects from defining elementary sends
 * @author Nicolas Sellenet
 */
class ObjectBalancer {
    /** @brief Vector of elementary sends to other processes */
    std::vector< VectorInt > _sendList;
    /** @brief Vector of number of elements to receive from others */
    VectorInt _recvSize;
    /** @brief Graph to browse */
    CommGraphPtr _graph;

    template < typename T >
    void balanceSimpleVectorOverProcesses( const T *, int, T * ) const;

  public:
    ObjectBalancer()
        : _sendList( std::vector< VectorInt >( getMPISize() ) ),
          _recvSize( VectorInt( getMPISize(), 0 ) ),
          _graph( new CommGraph() ){};

    /**
     * @brief Add an elementary send
     * @param rank Rank of receiver process
     * @param toSend vector of index to send
     */
    void addElementarySend( const int &rank, const VectorInt &toSend ) {
        _sendList[rank] = toSend;
        const auto result2 = std::min_element( toSend.begin(), toSend.end() );
        if ( *result2 < 0 )
            throw std::runtime_error( "Indexes of elements to send must be grower or equal to 0" );
        _graph->addCommunication( rank );
    };

    /** @brief Prepare communications (send and receive comm sizes) */
    void prepareCommunications();

    /** @brief Balance a object over processes by following elementary sends */
    // template< typename T >
    VectorReal balanceObjectOverProcesses( const VectorReal & ) const;
};

using ObjectBalancerPtr = std::shared_ptr< ObjectBalancer >;

template < typename T >
void ObjectBalancer::balanceSimpleVectorOverProcesses( const T *in, int sizeIn, T *out ) const {
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();
    VectorBool toKeep( sizeIn, true );
    for ( int iProc = 0; iProc < nbProcs; ++iProc ) {
        const auto curSendList = _sendList[iProc];
        const auto curSize = curSendList.size();
        for ( int iPos = 0; iPos < curSize; ++iPos ) {
            const auto vecPos = curSendList[iPos];
            if ( vecPos >= sizeIn )
                throw std::runtime_error( "Index to send grower to vector size" );
            if ( !toKeep[vecPos] )
                throw std::runtime_error( "Index " + std::to_string( vecPos ) +
                                          " in another comm (not allowed)" );
            toKeep[vecPos] = false;
        }
    }
    int cmpt = 0;
    for ( int iPos = 0; iPos < sizeIn; ++iPos ) {
        if ( toKeep[iPos] ) {
            out[cmpt] = in[iPos];
            ++cmpt;
        }
    }

    int tag = 0;
    for ( const auto proc : *_graph ) {
        if ( rank > proc ) {
            const auto curSendList = _sendList[proc];
            const auto curSendSize = curSendList.size();
            if ( curSendSize > 0 ) {
                VectorReal tmp( curSendSize, 0. );
                for ( int iPos = 0; iPos < curSendSize; ++iPos ) {
                    tmp[iPos] = in[curSendList[iPos]];
                }
                AsterMPI::send( tmp, proc, tag );
            }
            const auto curRecvSize = _recvSize[proc];
            if ( curRecvSize > 0 ) {
                VectorReal tmp( curRecvSize, 0. );
                AsterMPI::receive( tmp, proc, tag );
                for ( int curPos = 0; curPos < curRecvSize; ++curPos ) {
                    out[cmpt] = tmp[curPos];
                    ++cmpt;
                }
            }
        } else {
            const auto curRecvSize = _recvSize[proc];
            if ( curRecvSize > 0 ) {
                VectorReal tmp( curRecvSize, 0. );
                AsterMPI::receive( tmp, proc, tag );
                for ( int curPos = 0; curPos < curRecvSize; ++curPos ) {
                    out[cmpt] = tmp[curPos];
                    ++cmpt;
                }
            }
            const auto curSendList = _sendList[proc];
            const auto curSendSize = curSendList.size();
            if ( curSendSize > 0 ) {
                VectorReal tmp( curSendSize, 0. );
                for ( int iPos = 0; iPos < curSendSize; ++iPos ) {
                    tmp[iPos] = in[curSendList[iPos]];
                }
                AsterMPI::send( tmp, proc, tag );
            }
        }
        ++tag;
    }
};

#endif /* COMMGRAPH_H */
