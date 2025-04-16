#ifndef OBJECTBALANCER_H_
#define OBJECTBALANCER_H_

/**
 * @file ObjectBalancer.h
 * @brief Header of an object balancer
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
#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#include "aster_mpi.h"

#include "IOManager/MedVector.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxString.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "ParallelUtilities/CommGraph.h"
#include "ParallelUtilities/TemplateVectorTools.h"

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
    /** @brief Set of elements to keep on local process */
    std::set< int > _toKeep;
    /** @brief Set of elements to send */
    std::set< int > _toSend;
    /** @brief Set of elements to delete */
    std::set< int > _toDelete;
    /** @brief Size delta */
    int _sizeDelta;
    /** @brief Graph to browse */
    CommGraphPtr _graph;
    /** @brief Is objet ok to balance objects */
    bool _isOk;
    /** @brief To know if send definition has ended */
    bool _sendDefined;
    /** @brief Element renumbering */
    VectorLong _renumbering;

    template < typename T, int nbCmp = 1 >
    void balanceSimpleVectorOverProcesses( const T *, int, T * ) const;

  public:
    /**
     * @class DistributedMask
     * @brief Class used to apply a mask before sending and reverse after receiveing
     * @todo maybe a template ??
     * @author Nicolas Sellenet
     */
    class DistributedMask {
        const ObjectBalancer &_balancer;
        const VectorLong _vectMaskIn;
        const VectorLong _vectMaskOut;
        MapLong _mapMaskOut;

        MapLong buildMask() const {
            int cmpt = 1;
            MapLong mapMaskOut;
            for ( const auto globId : _vectMaskOut ) {
                mapMaskOut[globId] = cmpt;
                ++cmpt;
            }
            return mapMaskOut;
        }

      public:
        DistributedMask( const ObjectBalancer &balancer, const VectorLong &mask )
            : _balancer( balancer ),
              _vectMaskIn( mask ),
              _vectMaskOut( _balancer.balanceVectorOverProcesses( _vectMaskIn ) ),
              _mapMaskOut( buildMask() ) {};

        DistributedMask() = delete;

        const ASTERINTEGER &apply( const ASTERINTEGER &valueIn ) const {
            return _vectMaskIn[valueIn - 1];
        };

        const ASTERINTEGER &reverse( const ASTERINTEGER &valueIn ) const {
#ifdef ASTER_DEBUG_CXX
            if ( _mapMaskOut.find( valueIn ) == _mapMaskOut.end() ) {
                std::cout << "Num glob " << valueIn << " sans correspondance" << std::endl
                          << std::flush;
                throw std::runtime_error( "Error in DistributedMask " );
            }
#endif
            return _mapMaskOut.find( valueIn )->second;
        };

        const VectorLong &getBalancedMask() const { return _vectMaskOut; };
    };

    /**
     * @class DistributedMaskOut
     * @brief Class used to apply a mask only after receiveing
     * @todo maybe a template ??
     * @author Nicolas Sellenet
     */
    class DistributedMaskOut {
        const ObjectBalancer &_balancer;
        const VectorLong _vectMaskOut;
        MapLong _mapMaskOut;

        MapLong buildMask() const {
            int cmpt = 1;
            MapLong mapMaskOut;
            for ( const auto globId : _vectMaskOut ) {
                mapMaskOut[globId] = cmpt;
                ++cmpt;
            }
            return mapMaskOut;
        }

      public:
        DistributedMaskOut( const ObjectBalancer &balancer, const VectorLong &mask )
            : _balancer( balancer ),
              _vectMaskOut( _balancer.balanceVectorOverProcesses( mask ) ),
              _mapMaskOut( buildMask() ) {};

        DistributedMaskOut() = delete;

        ASTERINTEGER apply( const ASTERINTEGER &valueIn ) const { return valueIn; };

        const ASTERINTEGER &reverse( const ASTERINTEGER &valueIn ) const {
#ifdef ASTER_DEBUG_CXX
            if ( _mapMaskOut.find( valueIn ) == _mapMaskOut.end() ) {
                std::cout << "Num glob " << valueIn << " sans correspondance" << std::endl
                          << std::flush;
                throw std::runtime_error( "Error in DistributedMaskOut " );
            }
#endif
            return _mapMaskOut.find( valueIn )->second;
        };

        const VectorLong &getBalancedMask() const { return _vectMaskOut; };
    };

    struct DummyMask {
      public:
        ASTERINTEGER apply( const ASTERINTEGER &valueIn ) const { return valueIn; };

        ASTERINTEGER reverse( const ASTERINTEGER &valueIn ) const { return valueIn; };
    };

    struct DummyMaskDouble {
      public:
        double apply( const double &valueIn ) const { return valueIn; };

        double reverse( const double &valueIn ) const { return valueIn; };
    };

    template < typename Type >
    struct DummyMaskType {
      public:
        const Type apply( const Type &valueIn ) const { return valueIn; };

        const Type reverse( const Type &valueIn ) const { return valueIn; };
    };

  private:
    template < typename T, typename Mask = ObjectBalancer::DummyMask >
    void balanceObjectOverProcesses3( const T &, T &, const Mask &mask = Mask() ) const;

  public:
    ObjectBalancer()
        : _sendList( std::vector< VectorInt >( getMPISize() ) ),
          _recvSize( VectorInt( getMPISize(), 0 ) ),
          _sizeDelta( 0 ),
          _graph( std::make_shared< CommGraph >() ),
          _isOk( false ),
          _sendDefined( false ) {};

    /**
     * @brief Add an elementary send
     * @param rank Rank of receiver process
     * @param toSend vector of index to send
     */
    void addElementarySend( const int &rank, const VectorInt &toSend ) {
        if ( _sendDefined ) {
            throw std::runtime_error( "Definition of elementary sends already finished" );
        }
        if ( toSend.size() == 0 )
            return;
        _sendList[rank] = toSend;
        const auto result2 = std::min_element( toSend.begin(), toSend.end() );
        if ( *result2 < 0 )
            throw std::runtime_error( "Indexes of elements to send must be grower or equal to 0" );
        _graph->addCommunication( rank );
        for ( const auto &id : toSend ) {
            this->_toSend.insert( id );
        }
    };

    /**
     * @brief End elementary send definition
     */
    void endElementarySendDefinition() { _sendDefined = true; };

    /**
     * @brief Set list of elements to keep on local process (even if some are to send)
     * @param toKeep vector of index to keep
     */
    void setElementsToKeep( const VectorInt &toKeep ) {
        for ( const auto &id : toKeep ) {
            this->_toKeep.insert( id );
        }
    };

    /**
     * @brief Set list of elements to delete on local process
     * @param toDelete vector of index to delete
     */
    void setElementsToDelete( const VectorInt &toDelete ) {
        for ( const auto &id : toDelete ) {
            if ( _toKeep.find( id ) == _toKeep.end() )
                this->_toDelete.insert( id );
        }
    };

    /** @brief Prepare communications (send and receive comm sizes) */
    void prepareCommunications();

    /** @brief Function yo know if ObjectBalancer is useable (after prepareCommunications) */
    bool isUseable() const { return _isOk; };

    /** @brief Balance a object over processes by following elementary sends */
    template < typename T >
    T balanceVectorOverProcesses( const T & ) const;

    /** @brief Balance a object over processes by following elementary sends */
    template < typename T, int nbCmp = 1 >
    void balanceObjectOverProcesses( const T &, T & ) const;

    /** @brief Balance MeshCoordinatesField over processes by following elementary sends */
    void balanceObjectOverProcesses( const MeshCoordinatesFieldPtr &,
                                     MeshCoordinatesFieldPtr & ) const;

    /** @brief Balance JeveuxCollection over processes by following elementary sends */
    template < typename T, typename Mask = ObjectBalancer::DummyMask >
    void balanceObjectOverProcesses2( const JeveuxCollection< T, ASTERINTEGER, Contiguous > &,
                                      JeveuxCollection< T, ASTERINTEGER, Contiguous > &,
                                      const Mask &mask = Mask() ) const;
    template < typename T, typename Mask = ObjectBalancer::DummyMask >
    void balanceObjectOverProcesses2( const std::vector< std::vector< T > > &,
                                      std::vector< std::vector< T > > &,
                                      const Mask &mask = Mask() ) const;

    template < typename TypeName = double >
    std::shared_ptr< MedVector< TypeName > > balanceMedVectorOverProcessesWithRenumbering(
        const std::shared_ptr< MedVector< TypeName > > & ) const;

    template < typename TypeName = double >
    void
    balanceArrayOverProcessesWithRenumbering( const std::shared_ptr< ArrayWrapper< TypeName > > &,
                                              std::shared_ptr< ArrayWrapper< TypeName > > & ) const;

    VectorLong getRenumbering() const { return _renumbering; };

    void setRenumbering( const VectorLong &renumbering ) { _renumbering = renumbering; };
};

using ObjectBalancerPtr = std::shared_ptr< ObjectBalancer >;

template < typename T, int nbCmp >
void ObjectBalancer::balanceSimpleVectorOverProcesses( const T *in, int sizeIn, T *out ) const {
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();
    VectorBool toKeep( sizeIn, true );
    const auto sizeToKeep = _toKeep.size();
    const auto toKeepEnd = _toKeep.end();
    // Find elements to keep
    for ( int iProc = 0; iProc < nbProcs; ++iProc ) {
        const auto curSendList = _sendList[iProc];
        const auto curSize = curSendList.size();
        for ( int iPos = 0; iPos < curSize; ++iPos ) {
            for ( int iCmp = 0; iCmp < nbCmp; ++iCmp ) {
                const auto vecPos = nbCmp * curSendList[iPos] + iCmp;
                if ( vecPos >= sizeIn ) {
                    throw std::runtime_error( "Index to send grower than vector size " );
                }
                if ( sizeToKeep == 0 ) {
                    toKeep[vecPos] = false;
                } else {
                    if ( _toKeep.find( curSendList[iPos] ) == toKeepEnd )
                        toKeep[vecPos] = false;
                }
            }
        }
    }
    for ( const auto &id : _toDelete ) {
        for ( int iCmp = 0; iCmp < nbCmp; ++iCmp ) {
            const auto vecPos = nbCmp * id + iCmp;
            toKeep[vecPos] = false;
        }
    }
    // Copy of elements to keep in output vector
    int cmpt = 0;
    const auto curSize = sizeIn / nbCmp;
    for ( int iPos = 0; iPos < curSize; ++iPos ) {
        for ( int iCmp = 0; iCmp < nbCmp; ++iCmp ) {
            if ( toKeep[iPos * nbCmp + iCmp] ) {
                out[cmpt] = in[iPos * nbCmp + iCmp];
                ++cmpt;
            }
        }
    }

    // Loop over graph to communicate elements to send and receive
    // And then copy received values in output vector
    for ( const auto [tag, proc] : *_graph ) {
        if ( proc == -1 )
            continue;
        if ( rank > proc ) {
            const auto curSendList = _sendList[proc];
            const auto curSendSize = curSendList.size();
            if ( curSendSize > 0 ) {
                std::vector< T > tmp( curSendSize * nbCmp, 0. );
                for ( int iPos = 0; iPos < curSendSize; ++iPos ) {
                    for ( int iCmp = 0; iCmp < nbCmp; ++iCmp ) {
                        tmp[iPos * nbCmp + iCmp] = in[curSendList[iPos] * nbCmp + iCmp];
                    }
                }
                AsterMPI::send( tmp, proc, tag );
            }
            const auto curRecvSize = _recvSize[proc];
            if ( curRecvSize > 0 ) {
                std::vector< T > tmp( curRecvSize * nbCmp, 0. );
                AsterMPI::receive( tmp, proc, tag );
                for ( int curPos = 0; curPos < curRecvSize; ++curPos ) {
                    for ( int iCmp = 0; iCmp < nbCmp; ++iCmp ) {
                        out[cmpt] = tmp[curPos * nbCmp + iCmp];
                        ++cmpt;
                    }
                }
            }
        } else {
            const auto curRecvSize = _recvSize[proc];
            if ( curRecvSize > 0 ) {
                std::vector< T > tmp( curRecvSize * nbCmp, 0. );
                AsterMPI::receive( tmp, proc, tag );
                for ( int curPos = 0; curPos < curRecvSize; ++curPos ) {
                    for ( int iCmp = 0; iCmp < nbCmp; ++iCmp ) {
                        out[cmpt] = tmp[curPos * nbCmp + iCmp];
                        ++cmpt;
                    }
                }
            }
            const auto curSendList = _sendList[proc];
            const auto curSendSize = curSendList.size();
            if ( curSendSize > 0 ) {
                std::vector< T > tmp( curSendSize * nbCmp, 0. );
                for ( int iPos = 0; iPos < curSendSize; ++iPos ) {
                    for ( int iCmp = 0; iCmp < nbCmp; ++iCmp ) {
                        tmp[iPos * nbCmp + iCmp] = in[curSendList[iPos] * nbCmp + iCmp];
                    }
                }
                AsterMPI::send( tmp, proc, tag );
            }
        }
    }
};

template < typename T >
T ObjectBalancer::balanceVectorOverProcesses( const T &obj ) const {
    if ( !_isOk )
        throw std::runtime_error( "ObjectBalancer not prepared" );
    const auto vecSize = obj.size();
    T toReturn( vecSize + _sizeDelta, 0 );
    balanceSimpleVectorOverProcesses< typename T::value_type >( &obj[0], vecSize, &toReturn[0] );
    return toReturn;
};

template < typename T, int nbCmp >
void ObjectBalancer::balanceObjectOverProcesses( const T &in, T &out ) const {
    if ( !_isOk )
        throw std::runtime_error( "ObjectBalancer not prepared" );
    const auto vecSize = getSize< typename ValueType< T >::value_type >( (const T &)in );

    resize< typename ValueType< T >::value_type >( out, vecSize + _sizeDelta );
    balanceSimpleVectorOverProcesses< typename ValueType< T >::value_type, nbCmp >(
        getAddress( in ), vecSize, getAddress( out ) );
};

template < typename T, typename Mask >
void ObjectBalancer::balanceObjectOverProcesses2(
    const JeveuxCollection< T, ASTERINTEGER, Contiguous > &in,
    JeveuxCollection< T, ASTERINTEGER, Contiguous > &out, const Mask &mask ) const {
    balanceObjectOverProcesses3( *in, *out, mask );
};

template < typename T, typename Mask >
void ObjectBalancer::balanceObjectOverProcesses2( const std::vector< std::vector< T > > &in,
                                                  std::vector< std::vector< T > > &out,
                                                  const Mask &mask ) const {
    balanceObjectOverProcesses3( in, out, mask );
};

template < typename T, typename Mask >
void ObjectBalancer::balanceObjectOverProcesses3( const T &in, T &out, const Mask &mask ) const {
    if ( !_isOk )
        throw std::runtime_error( "ObjectBalancer not prepared" );
    typedef typename ObjectTemplateType< T >::value_type value_type;
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();
    const auto sizeIn = in.size();
    // if( sizeIn != 0 ) {
    //     if( !in->isContiguous() )
    //         throw std::runtime_error( "Only allowed for contiguous collection" );
    // }
    const auto totalSize = getTotalSize( in );
    VectorBool toKeep( sizeIn, true );
    const auto sizeToKeep = _toKeep.size();
    std::vector< VectorInt > occSize( nbProcs );
    VectorInt sizeToReceive( nbProcs, 0 );
    VectorInt sizeToSend( nbProcs, 0 );
    int toRemoveCum = 0, toReceiveCum = 0;
    constexpr int start = StartPosition< T >::value;
    for ( const auto [tag, proc] : *_graph ) {
        if ( proc == -1 )
            continue;
        const auto curSendList = _sendList[proc];
        const auto curSize = curSendList.size();
        int toRemove = 0, toSend = 0;
        for ( int iPos = 0; iPos < curSize; ++iPos ) {
            const auto vecPos = curSendList[iPos];
            const auto &occ = in[vecPos + start];
            const auto curSizeOC = getSize( occ );
            occSize[proc].push_back( curSizeOC );
            toSend += curSizeOC;
            if ( vecPos >= sizeIn ) {
                throw std::runtime_error( "Index to send grower than vector size" );
            }
            if ( sizeToKeep == 0 ) {
                if ( toKeep[vecPos] )
                    toRemove += curSizeOC;
                toKeep[vecPos] = false;
            } else {
                if ( _toKeep.find( curSendList[iPos] ) == _toKeep.end() ) {
                    if ( toKeep[vecPos] )
                        toRemove += curSizeOC;
                    toKeep[vecPos] = false;
                }
            }
        }
        toRemoveCum += toRemove;
        VectorInt tmp( 1, 0. );
        sizeToSend[proc] = toSend;
        if ( rank > proc ) {
            tmp[0] = toSend;
            AsterMPI::send( tmp, proc, tag );
            AsterMPI::receive( tmp, proc, tag );
            toReceiveCum += tmp[0];
            sizeToReceive[proc] = tmp[0];
        } else {
            AsterMPI::receive( tmp, proc, tag );
            toReceiveCum += tmp[0];
            sizeToReceive[proc] = tmp[0];
            tmp[0] = toSend;
            AsterMPI::send( tmp, proc, tag );
        }
    }
    for ( const auto &id : _toDelete ) {
        toKeep[id] = false;
    }

    allocate( out, sizeIn + _sizeDelta, totalSize - toRemoveCum + toReceiveCum );
    int cmpt = start;
    const auto curSize = sizeIn;
    for ( int iPos = 0; iPos < curSize; ++iPos ) {
        if ( toKeep[iPos] ) {
            auto toCopy = in[iPos + start];
            update( toCopy );
            const auto &sizetoCopy = getSize( toCopy );
            allocateOccurence( out, cmpt, sizetoCopy );
            auto &newObj = out[cmpt];
            for ( int curPos2 = 0; curPos2 < sizetoCopy; ++curPos2 ) {
                newObj[curPos2] = mask.reverse( mask.apply( toCopy[curPos2] ) );
            }
            ++cmpt;
        }
    }

    for ( const auto [tag, proc] : *_graph ) {
        if ( proc == -1 )
            continue;
        if ( rank > proc ) {
            const auto curSendList = _sendList[proc];
            const auto curSendSize = curSendList.size();
            if ( curSendSize > 0 ) {
                AsterMPI::send( occSize[proc], proc, tag );
                std::vector< value_type > tmp( sizeToSend[proc], 0. );
                int cmpt2 = 0;
                for ( int iPos = 0; iPos < curSendSize; ++iPos ) {
                    auto toCopy = in[curSendList[iPos] + start];
                    update( toCopy );
                    const auto &sizetoCopy = occSize[proc][iPos];
                    for ( int curPos2 = 0; curPos2 < sizetoCopy; ++curPos2 ) {
                        tmp[cmpt2] = mask.apply( toCopy[curPos2] );
                        ++cmpt2;
                    }
                }
                AsterMPI::send( tmp, proc, tag );
            }
            const auto curRecvSize = _recvSize[proc];
            if ( curRecvSize > 0 ) {
                VectorInt tmp( curRecvSize, 0. );
                AsterMPI::receive( tmp, proc, tag );
                std::vector< value_type > tmp2( sizeToReceive[proc], 0. );
                AsterMPI::receive( tmp2, proc, tag );
                int cmpt2 = 0;
                for ( int curPos = 0; curPos < curRecvSize; ++curPos ) {
                    const auto &sizetoCopy = tmp[curPos];
                    allocateOccurence( out, cmpt, sizetoCopy );
                    auto &newObj = out[cmpt];
                    for ( int curPos2 = 0; curPos2 < sizetoCopy; ++curPos2 ) {
                        newObj[curPos2] = mask.reverse( tmp2[cmpt2] );
                        ++cmpt2;
                    }
                    ++cmpt;
                }
            }
        } else {
            const auto curRecvSize = _recvSize[proc];
            if ( curRecvSize > 0 ) {
                VectorInt tmp( curRecvSize, 0. );
                AsterMPI::receive( tmp, proc, tag );
                std::vector< value_type > tmp2( sizeToReceive[proc], 0. );
                AsterMPI::receive( tmp2, proc, tag );
                int cmpt2 = 0;
                for ( int curPos = 0; curPos < curRecvSize; ++curPos ) {
                    const auto &sizetoCopy = tmp[curPos];
                    allocateOccurence( out, cmpt, sizetoCopy );
                    auto &newObj = out[cmpt];
                    for ( int curPos2 = 0; curPos2 < sizetoCopy; ++curPos2 ) {
                        newObj[curPos2] = mask.reverse( tmp2[cmpt2] );
                        ++cmpt2;
                    }
                    ++cmpt;
                }
            }
            const auto curSendList = _sendList[proc];
            const auto curSendSize = curSendList.size();
            if ( curSendSize > 0 ) {
                AsterMPI::send( occSize[proc], proc, tag );
                std::vector< value_type > tmp( sizeToSend[proc], 0. );
                int cmpt2 = 0;
                for ( int iPos = 0; iPos < curSendSize; ++iPos ) {
                    auto toCopy = in[curSendList[iPos] + start];
                    update( toCopy );
                    const auto &sizetoCopy = occSize[proc][iPos];
                    for ( int curPos2 = 0; curPos2 < sizetoCopy; ++curPos2 ) {
                        tmp[cmpt2] = mask.apply( toCopy[curPos2] );
                        ++cmpt2;
                    }
                }
                AsterMPI::send( tmp, proc, tag );
            }
        }
    }
};

template < typename TypeName >
std::shared_ptr< MedVector< TypeName > >
ObjectBalancer::balanceMedVectorOverProcessesWithRenumbering(
    const std::shared_ptr< MedVector< TypeName > > &vecIn ) const {
    std::shared_ptr< MedVector< TypeName > > vecOut( new MedVector< TypeName >() );
    balanceObjectOverProcesses3( *vecIn, *vecOut, DummyMaskDouble() );
    vecOut->setComponentNumber( vecIn->getComponentNumber() );
    vecOut->setComponentName( vecIn->getComponentName() );
    if ( _renumbering.size() == 0 ) {
        return vecOut;
    }
    const auto size = vecOut->size();
    std::shared_ptr< MedVector< TypeName > > vecOut2( new MedVector< TypeName >() );
    vecOut2->setComponentNumber( vecOut->getComponentNumber() );
    vecOut2->setComponentName( vecOut->getComponentName() );
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

template < typename TypeName >
void ObjectBalancer::balanceArrayOverProcessesWithRenumbering(
    const std::shared_ptr< ArrayWrapper< TypeName > > &vecIn,
    std::shared_ptr< ArrayWrapper< TypeName > > &vecOut ) const {
    typedef typename ArrayWrapper< TypeName >::value_type ValueType;
    // typedef typename std::shared_ptr< TypeName > TypeNamePtr;
    // std::shared_ptr< ArrayWrapper< TypeName > > vecOut( new ArrayWrapper< TypeName >() );
    balanceObjectOverProcesses3( *vecIn, *vecOut, DummyMaskType< ValueType >() );
    vecOut->setComponentNumber( vecIn->getComponentNumber() );
    vecOut->setComponentName( vecIn->getComponentName() );
    if ( _renumbering.size() != 0 ) {
        const auto size = vecOut->size();
        std::vector< ValueType > tmpVec;
        ArrayWrapper< std::vector< ValueType > > tmpArray =
            ArrayWrapper< std::vector< ValueType > >( tmpVec, vecOut->getComponentNumber() );
        tmpArray.setComponentNumber( vecOut->getComponentNumber() );
        tmpArray.setComponentName( vecOut->getComponentName() );
        tmpArray.setSize( size );
        if ( _renumbering.size() != size )
            throw std::runtime_error( "Sizes not matching" );
        for ( int i = 0; i < size; ++i ) {
            const auto newId = _renumbering[i] - 1;
            tmpArray.setElement( newId, vecOut->getElement( i ) );
        }
        tmpArray.endDefinition();
        for ( int i = 0; i < size; ++i ) {
            const auto newId = _renumbering[i] - 1;
            const auto nbCmp = vecOut->getElement( i );
            const auto &elInR = ( *vecOut )[i];
            auto &elOutR = tmpArray[newId];
            for ( int j = 0; j < nbCmp; ++j ) {
                elOutR[j] = elInR[j];
            }
        }
        vecOut->deallocate();
        vecOut->setComponentNumber( tmpArray.getComponentNumber() );
        vecOut->setComponentName( tmpArray.getComponentName() );
        vecOut->setSize( tmpArray.size() );
        for ( int i = 0; i < size; ++i ) {
            const auto newId = i;
            vecOut->setElement( newId, tmpArray.getElement( i ) );
        }
        vecOut->endDefinition();
        for ( int i = 0; i < size; ++i ) {
            const auto newId = i;
            const auto nbCmp = tmpArray.getElement( i );
            const auto &elInR = tmpArray[i];
            auto &elOutR = ( *vecOut )[newId];
            for ( int j = 0; j < nbCmp; ++j ) {
                elOutR[j] = elInR[j];
            }
        }
    }
};
#endif /* ASTER_HAVE_MPI */

#endif /* OBJECTBALANCER_H_ */
