/**
 * @file IncompleteMesh.cxx
 * @brief Implementation de IncompleteMesh
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

#include "Meshes/IncompleteMesh.h"

void IncompleteMesh::addFamily( int id, VectorString groups ) {
    if ( id > 0 ) {
        const auto size = _nodeFamGroups.size();
        if ( ( id - 1 ) >= size ) {
            for ( int i = size; i < id; ++i )
                _nodeFamGroups.push_back( VectorString() );
        }
        _nodeFamGroups[id - 1] = groups;
    } else if ( id < 0 ) {
        auto id2 = -id;
        const auto size = _cellFamGroups.size();
        if ( ( id2 - 1 ) >= size ) {
            for ( int i = size; i < id2; ++i )
                _cellFamGroups.push_back( VectorString() );
        }
        _cellFamGroups[id2 - 1] = groups;
    }
};

bool IncompleteMesh::debugCheckFromBaseMesh( BaseMeshPtr &compMesh ) const {
    bool toReturn = true;
    const auto &compCoords = *( compMesh->getCoordinates() );
    compCoords.updateValuePointers();
    const auto &coords = *_coordinates;
    coords.updateValuePointers();
    const auto nbNodes = coords.size() / 3;
    const auto firstId = _nodeRange[0];
    for ( int i = 0; i < nbNodes; ++i ) {
        const auto curNode = coords[i];
        const auto compCurNode = compCoords[firstId + i];
        if ( curNode[0] != compCurNode[0] ) {
#ifdef ASTER_DEBUG_CXX
            std::cout << "Diff X " << curNode[0] << " " << compCurNode[0] << std::endl;
#endif
            toReturn = false;
        }
        if ( curNode[1] != compCurNode[1] ) {
#ifdef ASTER_DEBUG_CXX
            std::cout << "Diff Y " << curNode[1] << " " << compCurNode[1] << std::endl;
#endif
            toReturn = false;
        }
        if ( curNode[2] != compCurNode[2] ) {
#ifdef ASTER_DEBUG_CXX
            std::cout << "Diff Z " << curNode[2] << " " << compCurNode[2] << std::endl;
#endif
            toReturn = false;
        }
    }

    auto &compConnex = *( compMesh->getConnectivity() );
    compConnex.build();
    auto &connex = *( getConnectivity() );
    connex.build();
    int cumElem = 0;
    for ( const auto &range : _cellRange ) {
        const auto &size = range[1] - range[0];
        for ( int localElemId = 0; localElemId < size; ++localElemId ) {
            auto &curConnex = connex[localElemId + cumElem + 1];
            curConnex->updateValuePointer();
            auto &curCheckConnex = compConnex[range[0] + localElemId + 1];
            curCheckConnex->updateValuePointer();
            const auto nbNode1 = curConnex->size();
            const auto nbNode2 = curCheckConnex->size();
            if ( nbNode1 != nbNode2 ) {
                toReturn = false;
                break;
            }
            for ( int j = 0; j < nbNode1; ++j ) {
                if ( ( *curConnex )[j] != ( *curCheckConnex )[j] ) {
                    toReturn = false;
                    break;
                }
            }
            if ( !toReturn )
                break;
        }
        cumElem += size;
        if ( !toReturn ) {
#ifdef ASTER_DEBUG_CXX
            std::cout << "Differences between connectivities" << std::endl;
#endif
            break;
        }
    }
    return toReturn;
};

ASTERINTEGER IncompleteMesh::getDimension() const {
    if ( isEmpty() )
        return 0;
    if ( !_dimensionInformations.exists() )
        return 0;

    _dimensionInformations->updateValuePointer();
    const auto dimGeom = ( *_dimensionInformations )[5];
    return dimGeom;
}

void IncompleteMesh::setCellFamily( const VectorLong &cf ) { _cellFamily = cf; };

VectorLong IncompleteMesh::getNodesFromGroup( const std::string &grpName ) const {
    VectorLong nodeIds;
    VectorLong foundFamId;
    int count = 1;
    for ( const auto &grpNames : _nodeFamGroups ) {
        for ( const auto &curName : grpNames ) {
            if ( curName == grpName ) {
                foundFamId.push_back( count );
            }
        }
        ++count;
    }
    for ( const auto &famIdToSearch : foundFamId ) {
        count = 0;
        for ( const auto &nodeFamId : _nodeFamily ) {
            if ( famIdToSearch == nodeFamId ) {
                nodeIds.push_back( count );
            }
            ++count;
        }
    }
    return nodeIds;
};

#endif /* ASTER_HAVE_MPI */
