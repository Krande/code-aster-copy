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

#include "Meshes/ParallelMesh.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"
#include "aster_fort_mesh.h"
#include "aster_fort_utils.h"

bool ParallelMesh::readPartitionedMedFile( const std::string &fileName ) {
    const bool ret = BaseMesh::readMedFile( fileName );

    CALLO_LRMJOI_WRAP( getName(), fileName );

    AS_ASSERT( updateGlobalGroupOfNodes() );
    AS_ASSERT( updateGlobalGroupOfCells() );

    return ret;
};

bool ParallelMesh::updateGlobalGroupOfNodes( void ) {

    _groupsOfNodes->build();
    auto gONNames = _groupsOfNodes->getObjectNames();
    std::vector< JeveuxChar32 > allgONNames;
    AsterMPI::all_gather( gONNames, allgONNames );

    for ( auto &nameOfGrp : allgONNames )
        _setOfAllGON.insert( trim( nameOfGrp.toString() ) );

    if ( _globalGroupOfNodes->exists() )
        _globalGroupOfNodes->deallocate();

    _globalGroupOfNodes->allocate( _setOfAllGON.size() );
    ASTERINTEGER num = 0;
    for ( auto &nameOfGrp : _setOfAllGON ) {
        _globalGroupOfNodes->add( num, nameOfGrp );
        ++num;
    }

    return true;
};

bool ParallelMesh::updateGlobalGroupOfCells( void ) {

    _groupsOfCells->build();
    auto gOENames = _groupsOfCells->getObjectNames();
    std::vector< JeveuxChar32 > allgOENames;
    AsterMPI::all_gather( gOENames, allgOENames );

    for ( auto &nameOfGrp : allgOENames )
        _setOfAllGOE.insert( trim( nameOfGrp.toString() ) );

    if ( _globalGroupOfCells->exists() )
        _globalGroupOfCells->deallocate();

    _globalGroupOfCells->allocate( _setOfAllGOE.size() );
    ASTERINTEGER num = 0;
    for ( auto &nameOfGrp : _setOfAllGOE ) {
        _globalGroupOfCells->add( num, nameOfGrp );
        ++num;
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
        return irange( long( 1 ), long( getNumberOfCells() ) );
    } else if ( !hasGroupOfCells( name, true ) ) {
        return VectorLong();
    }

    return _groupsOfCells->getObjectFromName( name ).toVector();
}

VectorLong ParallelMesh::getNodes( const std::string name, const bool localNumbering ) const {
    CALL_JEMARQ();
    VectorLong listOfNodes;
    if ( name.empty() ) {
        listOfNodes = irange( long( 1 ), long( getNumberOfNodes() ) );
    } else if ( !hasGroupOfNodes( name, true ) ) {
        return VectorLong();
    } else {
        listOfNodes = _groupsOfNodes->getObjectFromName( name ).toVector();
    }

    if ( localNumbering )
        return listOfNodes;

    VectorLong newNumbering;
    newNumbering.reserve( listOfNodes.size() );
    _globalNumbering->updateValuePointer();

    for ( auto &nodeId : listOfNodes )
        newNumbering.push_back( ( *_globalNumbering )[nodeId - 1] );
    CALL_JEDEMA();

    return newNumbering;
}

VectorLong ParallelMesh::getNodes( const std::string name, const bool localNumbering,
                                   const bool same_rank ) const {

    CALL_JEMARQ();
    VectorLong listOfNodes;
    if ( name.empty() ) {
        listOfNodes = irange( long( 1 ), long( getNumberOfNodes() ) );
    } else if ( !hasGroupOfNodes( name, true ) ) {
        return VectorLong();
    } else {
        listOfNodes = _groupsOfNodes->getObjectFromName( name ).toVector();
    }

    const int rank = getMPIRank();
    VectorLong newRank;
    newRank.reserve( listOfNodes.size() );

    int size = 0;
    _outerNodes->updateValuePointer();
    for ( int nodeId : listOfNodes ) {
        if ( same_rank ) {
            if ( rank == ( *_outerNodes )[nodeId - 1] ) {
                newRank.push_back( nodeId );
                size++;
            }
        } else {
            if ( rank != ( *_outerNodes )[nodeId - 1] ) {
                newRank.push_back( nodeId );
                size++;
            }
        }
    }
    newRank.resize( size );

    if ( localNumbering )
        return newRank;

    VectorLong newNumbering;
    newNumbering.reserve( newRank.size() );
    _globalNumbering->updateValuePointer();
    for ( auto &nodeId : newRank )
        newNumbering.push_back( ( *_globalNumbering )[nodeId - 1] );

    CALL_JEDEMA();

    return newNumbering;
};

VectorLong ParallelMesh::getNodesFromCells( const std::string name, const bool localNumbering,
                                            const bool same_rank ) const {
    CALL_JEMARQ();
    const auto cellsId = getCells( name );

    if ( cellsId.empty() )
        return VectorLong();

    const auto &connecExp = getConnectivityExplorer();

    SetLong nodes;

    for ( auto &cellId : cellsId ) {
        const auto cell = connecExp[cellId];
        for ( auto &node : cell )
            nodes.insert( node );
    }

    if ( same_rank ) {
        _outerNodes->updateValuePointer();
        const int rank = getMPIRank();

        for ( auto &node : nodes ) {
            if ( rank != ( *_outerNodes )[node - 1] )
                nodes.erase( node );
        }
    }

    if ( !localNumbering ) {
        VectorLong v_nodes;
        v_nodes.reserve( nodes.size() );

        _globalNumbering->updateValuePointer();
        for ( auto &node : nodes )
            v_nodes.push_back( ( *_globalNumbering )[node - 1] );

        return v_nodes;
    }

    CALL_JEDEMA();

    return VectorLong( nodes.begin(), nodes.end() );
};

bool ParallelMesh::build(){

    CALL_JEMARQ();
    _listOfOppositeDomain->updateValuePointer();
    const auto size = _listOfOppositeDomain->size();

    const std::string cadre("G");

    for(ASTERINTEGER i = 0; i < size; i++ ){
        auto domdis = (*_listOfOppositeDomain)[i];
        std::string chdomdis(8, ' ');

        CALLO_CODENT(&domdis, cadre, chdomdis);
        //std::cout << chdomdis << std::endl;
        JeveuxVectorLong jointE(getName() + ".E" + chdomdis);
        JeveuxVectorLong jointR(getName() + ".R" + chdomdis);

        _joints[domdis] = std::make_pair(jointE, jointR);
    }

    CALL_JEDEMA();

    return BaseMesh::build();
}

#endif /* ASTER_HAVE_MPI */
