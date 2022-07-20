/**
 * @file ParallelMesh.cxx
 * @brief Implementation de ParallelMesh
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */
#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#include "aster_fort_mesh.h"
#include "aster_fort_utils.h"

#include "Meshes/ParallelMesh.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

bool ParallelMesh::readPartitionedMedFile( const std::string &fileName ) {
    const bool ret = BaseMesh::readMedFile( fileName );

    CALLO_LRMJOI_WRAP( getName(), fileName );

    AS_ASSERT( updateGlobalGroupOfNodes() );
    AS_ASSERT( updateGlobalGroupOfCells() );
    AS_ASSERT( build() );

    return ret;
};

bool ParallelMesh::updateGlobalGroupOfNodes( void ) {

    _groupsOfNodes->build();
    auto gONNames = _groupsOfNodes->getObjectsNames();
    std::vector< JeveuxChar32 > allgONNames;
    AsterMPI::all_gather( gONNames, allgONNames );

    for ( auto &nameOfGrp : allgONNames )
        _setOfAllGON.insert( trim( nameOfGrp.toString() ) );

    if ( _globalGroupOfNodes->exists() )
        _globalGroupOfNodes->deallocate();

    if ( _setOfAllGON.size() > 0 ) {
        _globalGroupOfNodes->allocate( _setOfAllGON.size() );
        ASTERINTEGER num = 0;
        for ( auto &nameOfGrp : _setOfAllGON ) {
            _globalGroupOfNodes->add( num, nameOfGrp );
            ++num;
        }
    }

    return true;
};

bool ParallelMesh::updateGlobalGroupOfCells( void ) {

    _groupsOfCells->build();
    auto gOENames = _groupsOfCells->getObjectsNames();
    std::vector< JeveuxChar32 > allgOENames;
    AsterMPI::all_gather( gOENames, allgOENames );

    for ( auto &nameOfGrp : allgOENames )
        _setOfAllGOE.insert( trim( nameOfGrp.toString() ) );

    if ( _globalGroupOfCells->exists() )
        _globalGroupOfCells->deallocate();

    if ( _setOfAllGOE.size() > 0 ) {
        _globalGroupOfCells->allocate( _setOfAllGOE.size() );
        ASTERINTEGER num = 0;
        for ( auto &nameOfGrp : _setOfAllGOE ) {
            _globalGroupOfCells->add( num, nameOfGrp );
            ++num;
        }
    }

    return true;
};

bool ParallelMesh::hasGroupOfCells( const std::string &name, const bool local ) const {

    if ( local ) {
        if ( _groupsOfCells->size() < 0 && !_groupsOfCells->build() ) {
            return false;
        }

        return _groupsOfCells->existsObject( name );
    }

    SetOfStringCIter curIter = _setOfAllGOE.find( name );
    if ( curIter != _setOfAllGOE.end() )
        return true;
    return false;
}

bool ParallelMesh::hasGroupOfNodes( const std::string &name, const bool local ) const {
    if ( local ) {
        if ( _groupsOfNodes->size() < 0 && !_groupsOfNodes->build() ) {
            return false;
        }
        return _groupsOfNodes->existsObject( name );
    }

    SetOfStringCIter curIter = _setOfAllGON.find( name );
    if ( curIter != _setOfAllGON.end() )
        return true;
    return false;
}

VectorString ParallelMesh::getGroupsOfCells( const bool local ) const {

    if ( local ) {
        ASTERINTEGER size = _nameOfGrpCells->size();
        VectorString names;
        for ( int i = 0; i < size; i++ ) {
            names.push_back( trim( _nameOfGrpCells->getStringFromIndex( i + 1 ) ) );
        }
        return names;
    }

    return VectorString( _setOfAllGOE.begin(), _setOfAllGOE.end() );
}

VectorString ParallelMesh::getGroupsOfNodes( const bool local ) const {

    if ( local ) {
        ASTERINTEGER size = _nameOfGrpNodes->size();
        VectorString names;
        for ( int i = 0; i < size; i++ ) {
            names.push_back( trim( _nameOfGrpNodes->getStringFromIndex( i + 1 ) ) );
        }
        return names;
    }

    return VectorString( _setOfAllGON.begin(), _setOfAllGON.end() );
}

VectorLong ParallelMesh::getCells( const std::string name ) const {

    if ( name.empty() ) {
        return irange( long( 0 ), long( getNumberOfCells() - 1 ) );
    } else if ( !hasGroupOfCells( name, true ) ) {
        return VectorLong();
    }

    VectorLong cells = ( *_groupsOfCells )[name].toVector();
    for ( auto &cell : cells )
        cell -= 1;
    return cells;
}

VectorLong ParallelMesh::getNodes( const std::string name, const bool localNumbering,
                                   const ASTERINTEGER same_rank ) const {
    CALL_JEMARQ();
    VectorLong listOfNodes;
    if ( name.empty() ) {
        listOfNodes = irange( long( 0 ), long( getNumberOfNodes() - 1 ) );
    } else if ( !hasGroupOfNodes( name, true ) ) {
        CALL_JEDEMA();
        return VectorLong();
    } else {
        listOfNodes = ( *_groupsOfNodes )[name].toVector();
        for ( auto &node : listOfNodes )
            node -= 1;
    }

    const int rank = getMPIRank();
    VectorLong newRank;

    if ( same_rank == PythonBool::None ) {
        newRank = listOfNodes;
    } else {
        newRank.reserve( listOfNodes.size() );
        ASTERINTEGER size = 0;
        _outerNodes->updateValuePointer();
        for ( auto &nodeId : listOfNodes ) {
            if ( same_rank ) {
                if ( rank == ( *_outerNodes )[nodeId] ) {
                    newRank.push_back( nodeId );
                    size++;
                }
            } else {
                if ( rank != ( *_outerNodes )[nodeId] ) {
                    newRank.push_back( nodeId );
                    size++;
                }
            }
        }
        newRank.resize( size );
    }

    listOfNodes.clear();

    if ( localNumbering ) {
        CALL_JEDEMA();
        return newRank;
    }

    VectorLong newNumbering;
    newNumbering.reserve( newRank.size() );
    _globalNumbering->updateValuePointer();

    for ( auto &nodeId : newRank )
        newNumbering.push_back( ( *_globalNumbering )[nodeId] );
    CALL_JEDEMA();

    return newNumbering;
}

VectorLong ParallelMesh::getNodesFromCells( const VectorLong &cells, const bool localNumbering,
                                            const ASTERINTEGER same_rank ) const {
    if ( cells.empty() ) {
        return VectorLong();
    }

    CALL_JEMARQ();

    const auto &connecExp = getConnectivityExplorer();

    SetLong nodes;

    for ( auto &cellId : cells ) {
        const auto cell = connecExp[cellId];
        for ( auto &node : cell )
            nodes.insert( node - 1 );
    }

    if ( same_rank != PythonBool::None ) {
        _outerNodes->updateValuePointer();
        const int rank = getMPIRank();

        SetLong loopnodes = nodes;
        for ( auto &node : loopnodes ) {
            if ( same_rank ) {
                if ( rank != ( *_outerNodes )[node] )
                    nodes.erase( node );
            } else {
                if ( rank == ( *_outerNodes )[node] )
                    nodes.erase( node );
            }
        }
    }

    if ( !localNumbering ) {
        VectorLong v_nodes;
        v_nodes.reserve( nodes.size() );

        _globalNumbering->updateValuePointer();
        for ( auto &node : nodes )
            v_nodes.push_back( ( *_globalNumbering )[node] );

        CALL_JEDEMA();
        return v_nodes;
    }

    CALL_JEDEMA();

    return VectorLong( nodes.begin(), nodes.end() );
};

VectorLong ParallelMesh::getNodesFromCells( const std::string name, const bool localNumbering,
                                            const ASTERINTEGER same_rank ) const {
    return getNodesFromCells( getCells( name ), localNumbering, same_rank );
};

VectorLong ParallelMesh::getInnerCells() const {
    CALL_JEMARQ();
    VectorLong listOfCells = getCells();

    const int rank = getMPIRank();
    VectorLong newRank;

    newRank.reserve( listOfCells.size() );
    _outerCells->updateValuePointer();

    for ( auto &cell : listOfCells ) {
        if ( rank == ( *_outerCells )[cell] ) {
            newRank.push_back( cell );
        }
    }

    newRank.shrink_to_fit();
    CALL_JEDEMA();

    return newRank;
}

VectorLong ParallelMesh::getOuterCells() const {
    CALL_JEMARQ();
    VectorLong listOfCells = getCells();

    const int rank = getMPIRank();
    VectorLong newRank;

    newRank.reserve( listOfCells.size() );
    _outerCells->updateValuePointer();

    for ( auto &cell : listOfCells ) {
        if ( rank != ( *_outerCells )[cell] ) {
            newRank.push_back( cell );
        }
    }

    newRank.shrink_to_fit();
    CALL_JEDEMA();

    return newRank;
}

bool ParallelMesh::build() {

    const auto size = _listOfOppositeDomain->size();
    if ( _joints.empty() && size > 0 ) {
        CALL_JEMARQ();
        _listOfOppositeDomain->updateValuePointer();

        const std::string cadre( "G" );
        const std::string error( "F" );

        for ( ASTERINTEGER i = 0; i < size; i++ ) {
            auto domdis = ( *_listOfOppositeDomain )[i];
            std::string chdomdis( 4, ' ' );

            CALLO_CODLET_WRAP( &domdis, cadre, chdomdis, error );
            JeveuxVectorLong jointE( getName() + ".E" + chdomdis );
            JeveuxVectorLong jointR( getName() + ".R" + chdomdis );

            _joints[domdis] = std::make_pair( jointE, jointR );
        }

        CALL_JEDEMA();
    }

    return BaseMesh::build();
}

#endif /* ASTER_HAVE_MPI */
