/**
 * @file ConnectionMesh.cxx
 * @brief Implementation de
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

#include <algorithm>
#include <sstream>
#include <iomanip>

#include "aster_fort_mesh.h"
#include "Meshes/ConnectionMesh.h"
#include "ParallelUtilities/AsterMPI.h"

#ifdef ASTER_HAVE_MPI

/* Initial constructor for an object of ConnectionMesh. This class is only
   defined for the parallel use of a ParrallelMesh. The idea is to rebuild on
   each processor a temporary mesh associated with all the cells and nodes
   concerned by the groups given.
*/
ConnectionMesh::ConnectionMesh( const std::string &name,
                                          const ParallelMeshPtr &mesh,
                                          const VectorString &groupsOfNodes,
                                          const VectorString &groupsOfCells )
: BaseMesh( name, "MAILLAGE_PARTIEL" ),
  _pMesh( mesh ),
  _nodesLocalNumbering( getName() + ".NOLOCAL" ),
  _nodesGlobalNumbering( getName() + ".NOGLOBAL" ),
  _nodesOwner( getName() + ".NOPROPRIO" ),
  _cellsLocalNumbering( getName() + ".MALOCAL"),
  _cellsOwner( getName() + ".MAPROPRIO")
{

    /* Stop if no group is given */
    if( groupsOfNodes.empty() && groupsOfCells.empty() )
        throw std::runtime_error( "No groups" );

    /* ******************************************************************* */
    /* Part of the code dealing with the definition of the local variables */
    /* ******************************************************************* */

    /* Local variables gathering the basic MPI informations needed */
    /* Rank of the current MPI proc */
    const int rank = getMPIRank();

    /* Local variables needed to handle the input groups of cells */
     /* Make a copy of the input groups of cells */
    VectorString groupsOfCellsToFind;
    groupsOfCellsToFind.reserve( groupsOfCells.size() );
    /* Total number of mesh cells related to the proc */
    const int numberOfMeshCells = mesh->getNumberOfCells();
    /* Rank of the proc cells (to identify outer cells) */
    const JeveuxVectorLong rankOfCells = mesh->getCellsRank();
    /* Update the jeveux pointer to access the Fortran Jeveux pointer */
    rankOfCells->updateValuePointer();
    /* Get global cell connectivity to create cells group's later */
    const VectorLong globalCellIds = this->getCellsGlobalNumbering( rankOfCells );
    int numberOfCellsToSend = 0;
    /* Cells of the proc that are part of the ConnectionMesh in building */
    VectorLong cellsToSend;
    cellsToSend.reserve( numberOfMeshCells );
    /* Tags related to proc cells, used to build their future numbering */
    VectorBool boolCellsToSend( numberOfMeshCells, false );
    /* Outer cells whose future numbering needs checking in other procs */
    int numberOfCellsToCheck = 0;
    VectorLong cellsToCheck;
    cellsToCheck.reserve( numberOfMeshCells );

    std::map< std::string, VectorLong > groupsOfCellsToSend;
    std::map< std::string, VectorLong > groupsOfCellsGathered;

    VectorLong connectivitiesToSend;
    VectorLong connectivitiesGathered;

    /* Local variables needed to handle the input groups of nodes */
    /* Make a copy of the input groups of nodes */
    VectorString groupsOfNodesToFind;
    groupsOfNodesToFind.reserve( groupsOfNodes.size() );
    /* Total number of mesh nodes related to the proc */
    const int numberOfMeshNodes = mesh->getNumberOfNodes();
    /* Rank of the proc nodes (to identify outer nodes) */
    const JeveuxVectorLong rankOfNodes = mesh->getNodesRank();
    /* Update the jeveux pointer to access the Fortran Jeveux pointer */
    rankOfNodes->updateValuePointer();
    /* Global numbering of proc nodes (to link outer nodes on all procs) */
    const VectorLong globalNodeIds = mesh->getNodes( false );
    int numberOfNodesToSend = 0;
    /* The proc Nodes that are part of the ConnectionMesh we are building */
    VectorLong nodesToSend;
    nodesToSend.reserve( numberOfMeshNodes );
    VectorLong outerNodesToSend;
    outerNodesToSend.reserve( mesh->getOuterNodes().size() );

    /* Tags related to the proc nodes used to build future numbering */
    VectorBool boolNodesToSend( numberOfMeshNodes, false );
    /* Outer nodes whose future numbering needs checking in other procs */
    int numberOfNodesToCheck = 0;
    VectorLong nodesToCheck;
    nodesToCheck.reserve( numberOfMeshNodes );

    std::map< std::string, VectorLong > groupsOfNodesToSend;
    std::map< std::string, VectorLong > groupsOfNodesGathered;

    const auto &meshCoordinates = mesh->getCoordinates();
    meshCoordinates->updateValuePointers();

    /* Inverse connectivity nodes -> cell */
    const JeveuxCollectionLong connecInv = mesh->getInverseConnectivity();
    connecInv->build();

    VectorReal coordinatesToSend;
    VectorReal coordinatesGathered;

    VectorLong numNodesGathered;
    VectorLong numNodesToSend;

    std::map< long, long> numNodesGloLoc, numCellsGloLoc, renumNodes;

    int totalNumberOfNodes = 0, totalNumberOfCells = 0, numOwner = 0;
    long pos = 0;

    const auto& connecExp = mesh->getConnectivityExplorer();

    /* *********************************************************************** */
    /* Part of the code dealing with the treament of the input groups of cells */
    /* *********************************************************************** */

    /* Sort copy of input vector to be identical on all procs */
    for(const auto &nameOfTheGroup : groupsOfCells )
    {
        if( mesh->hasGroupOfCells( nameOfTheGroup, false ) )
            groupsOfCellsToFind.push_back( nameOfTheGroup );
    }
    std::sort( groupsOfCellsToFind.begin(), groupsOfCellsToFind.end() );

    if ( !groupsOfCellsToFind.empty() )
    {
        /* Loop over the groups of cells to add/tag those concerned by the proc,
        add/tag also the nodes concerned and get the ones needing later checks */
        for ( const auto &nameOfTheGroup : groupsOfCellsToFind )
        {
            if( mesh->hasGroupOfCells( nameOfTheGroup, true ) )
            {
                const auto &cellsToFind = mesh->getCells( nameOfTheGroup );
                const auto numberOfCellsToFind = cellsToFind.size();
                VectorLong cellsOfTheGroupToSend;
                cellsOfTheGroupToSend.reserve( numberOfCellsToFind );
                for ( int i = 0; i < numberOfCellsToFind; ++i )
                {
                    const auto cellId = cellsToFind[i] - 1;

                    if ( !boolCellsToSend[cellId] )
                    {
                        const auto cell = connecExp[cellsToFind[i]];
                        for (const auto vertex : cell )
                        {
                            const auto nodeId = vertex - 1;
                            if ( !boolNodesToSend[nodeId] )
                            {
                                /*  Split the proc nodes (toSend) from outer ones (toCheck) */
                                if ( ( *rankOfNodes )[nodeId] == rank )
                                {
                                    nodesToSend.push_back( vertex );
                                    ++numberOfNodesToSend;
                                }
                                else
                                {
                                    nodesToCheck.push_back( vertex );
                                    ++numberOfNodesToCheck;
                                }
                                boolNodesToSend[nodeId] = true;
                            }
                        }
                    }

                    /* Add the cell that we own */
                    if ( ( *rankOfCells )[cellId] == rank )
                    {
                        cellsOfTheGroupToSend.push_back( globalCellIds[cellId] );

                        if ( !boolCellsToSend[cellId] )
                        {
                            cellsToSend.push_back( cellsToFind[i] );
                            ++numberOfCellsToSend;
                        }
                    }
                    boolCellsToSend[cellId] = true;
                }
                groupsOfCellsToSend[nameOfTheGroup] = cellsOfTheGroupToSend;
            }
            groupsOfCellsGathered[nameOfTheGroup] = VectorLong();
        }
    }


    /* *********************************************************************** */
    /* Part of the code dealing with the treament of the input groups of nodes */
    /* *********************************************************************** */
    /* Sort copy of input vector to be identical on all procs */
    for(const auto &nameOfTheGroup : groupsOfNodes )
    {
        if( mesh->hasGroupOfNodes( nameOfTheGroup, false ) )
            groupsOfNodesToFind.push_back( nameOfTheGroup );
    }
    std::sort( groupsOfNodesToFind.begin(), groupsOfNodesToFind.end() );

    if ( !groupsOfNodesToFind.empty() )
    {
        /* Loop over the groups of nodes to tag/add those concerned by the proc */
        for ( const auto &nameOfTheGroup : groupsOfNodesToFind )
        {
            if ( mesh->hasGroupOfNodes( nameOfTheGroup, true ) )
            {
                const auto &nodesToFind = mesh->getNodes( nameOfTheGroup );
                const int numberOfNodesToFind = nodesToFind.size();
                VectorLong nodesOfTheGroupToSend;
                nodesOfTheGroupToSend.reserve(numberOfNodesToFind);
                for ( int i = 0; i < numberOfNodesToFind; ++i )
                {
                    const auto nodeId = nodesToFind[i] - 1;
                    if ( !boolNodesToSend[nodeId] )
                    {
                        if ( ( *rankOfNodes )[nodeId] == rank )
                        {
                            nodesToSend.push_back( nodesToFind[i] );
                            ++numberOfNodesToSend;
                        }
                        else
                        {
                            nodesToCheck.push_back( nodesToFind[i] );
                            ++numberOfNodesToCheck;
                        }
                        boolNodesToSend[nodeId] = true;
                    }

                    if ( ( *rankOfNodes )[nodeId] == rank )
                    {
                        nodesOfTheGroupToSend.push_back( globalNodeIds[nodeId] );
                    }
                }
                groupsOfNodesToSend[nameOfTheGroup] = nodesOfTheGroupToSend;
            }
            groupsOfNodesGathered[nameOfTheGroup] = VectorLong();
        }
    }


    /* *********************************************************************** */
    /* Part of the code dealing with the neighbourhood                         */
    /* *********************************************************************** */

    /* Get cells lying on nodes */
    for( int i = 0; i < numberOfNodesToSend; i++)
    {
        const auto nodeId = nodesToSend[i] - 1;
        const auto listCells = connecInv->getObject( nodeId + 1 ).toVector();

        for( const auto cell : listCells )
        {
            const auto cellId = cell - 1;
            if ( !boolCellsToSend[cellId] )
            {
                if ( ( *rankOfCells )[cellId] == rank )
                {
                    cellsToSend.push_back( cell );
                    ++numberOfCellsToSend;
                }
                else
                {
                    cellsToCheck.push_back( cell );
                    ++numberOfCellsToCheck;
                }
                boolCellsToSend[cellId] = true;
            }
        }
    }

    for( int i = 0; i < numberOfNodesToCheck; i++)
    {
        const auto nodeId = nodesToCheck[i] - 1;
        const auto listCells = connecInv->getObject( nodeId + 1 ).toVector();

        for( const auto cell : listCells )
        {
            const auto cellId = cell - 1;
            if ( !boolCellsToSend[cellId] )
            {
                if ( ( *rankOfCells )[cellId] == rank )
                {
                    cellsToSend.push_back( cell );
                    ++numberOfCellsToSend;
                }
                else
                {
                    cellsToCheck.push_back( cell );
                    ++numberOfCellsToCheck;
                }
                boolCellsToSend[cellId] = true;
            }
        }
    }
    nodesToCheck.clear();

    /* Get nodes lying on cell */
    for( int i = 0; i < numberOfCellsToSend; i++)
    {
        const auto cell = connecExp[cellsToSend[i]];
        for (const auto vertex : cell )
        {
            const auto nodeId = vertex - 1;
            if( !boolNodesToSend[nodeId] )
            {
                if ( ( *rankOfNodes )[nodeId] == rank )
                {
                    nodesToSend.push_back( vertex );
                    ++numberOfNodesToSend;
                }
                else
                {
                    outerNodesToSend.push_back(globalNodeIds[nodeId]);
                    outerNodesToSend.push_back((*rankOfNodes)[nodeId]);
                }
            }
            boolNodesToSend[nodeId] = true;
        }
    }

    for( int i = 0; i < numberOfCellsToCheck; i++)
    {
        const auto cell = connecExp[cellsToCheck[i]];
        for (const auto vertex : cell )
        {
            const auto nodeId = vertex - 1;
            if ( !boolNodesToSend[nodeId] && ( *rankOfNodes )[nodeId] == rank )
            {
                nodesToSend.push_back( vertex );
                ++numberOfNodesToSend;
            }
            boolNodesToSend[nodeId] = true;
        }
    }
    cellsToCheck.clear();

    /* Some nodes can be not marked if they are not owned by the current proc
    * We have to collect and mark them
    * */

    outerNodesToSend.shrink_to_fit();
    VectorLong outerNodesToGathered;
    AsterMPI::all_gather( outerNodesToSend, outerNodesToGathered );
    outerNodesToSend.clear();

    std::map< long, long > inverseGlobalNodeIds;
    pos = 0;
    for(auto globalId : globalNodeIds)
    {
        inverseGlobalNodeIds[globalId] = pos++;
    }

    const auto nbOuterNodes = outerNodesToGathered.size() / 2;
    for(int i = 0; i < nbOuterNodes; i++)
    {
        const auto globalNodeId = outerNodesToGathered[2*i];
        const auto ownerRank = outerNodesToGathered[2*i+1];
        if( rank == ownerRank )
        {
            const auto localNodeId = inverseGlobalNodeIds[globalNodeId];
            if( !boolNodesToSend[localNodeId] )
            {
                nodesToSend.push_back( localNodeId + 1 );
                ++numberOfNodesToSend;
                boolNodesToSend[localNodeId] = true;
            }
        }
    }

    outerNodesToGathered.clear();
    inverseGlobalNodeIds.clear();

    boolNodesToSend.clear();
    boolCellsToSend.clear();

    /* Now we have all cells and nodes to send */


    /* *********************************************************************** */
    /* Part of the code dealing with the building of the output ConnectionMesh */
    /* *********************************************************************** */

    /* Build the coordinates and numbering of the mesh nodes to send */
    coordinatesToSend.reserve( 3 * nodesToSend.size());
    numNodesToSend.reserve( 3 * nodesToSend.size());
    for ( const auto &node : nodesToSend )
    {
        const auto nodeId = node - 1;
        coordinatesToSend.push_back( ( *meshCoordinates )[nodeId * 3] );
        coordinatesToSend.push_back( ( *meshCoordinates )[nodeId * 3 + 1] );
        coordinatesToSend.push_back( ( *meshCoordinates )[nodeId * 3 + 2] );
        numNodesToSend.push_back( node );
        numNodesToSend.push_back( globalNodeIds[nodeId] );
        numNodesToSend.push_back( rank );
    }
    nodesToSend.clear();

    AsterMPI::all_gather( coordinatesToSend, coordinatesGathered );
    totalNumberOfNodes = coordinatesGathered.size() / 3;
    coordinatesToSend.clear();

    AsterMPI::all_gather( numNodesToSend, numNodesGathered );
    numNodesToSend.clear();


    /* Build the connectivities of the mesh cells and their types to send */
    /* For each cell, we save the type, the global index, the number of nodes,
        the list of nodes in global numbering */
    connectivitiesToSend.reserve( cellsToSend.size() * ( 1 + 1 + 1 + 27 ));
    for ( const auto& cellId : cellsToSend )
    {
        const auto cell = connecExp[cellId];
        connectivitiesToSend.push_back( cell.getType() );
        connectivitiesToSend.push_back( cellId );
        connectivitiesToSend.push_back( globalCellIds[ cellId - 1 ] );
        connectivitiesToSend.push_back( rank );
        connectivitiesToSend.push_back( cell.getNumberOfNodes() );
        for (const auto vertex : cell )
        {
            const auto nodeId = vertex - 1;
            connectivitiesToSend.push_back( globalNodeIds[nodeId] );
        }
    }

    /* Gather the types and connectivities of cells */
    AsterMPI::all_reduce(int(cellsToSend.size()), totalNumberOfCells, MPI_SUM);
    cellsToSend.clear();

    AsterMPI::all_gather( connectivitiesToSend, connectivitiesGathered );
    connectivitiesToSend.clear();

    /* Gather the groups of nodes */
    for ( const auto &nameOfTheGroup : groupsOfNodesToFind )
    {
        VectorLong &nodesOfTheGroupToSend = groupsOfNodesToSend[nameOfTheGroup];
        VectorLong &nodesOfTheGroupGathered = groupsOfNodesGathered[nameOfTheGroup];

        AsterMPI::all_gather( nodesOfTheGroupToSend, nodesOfTheGroupGathered );
        nodesOfTheGroupToSend.clear();
    }

    /* Gather the groups of cells */
    for ( const auto &nameOfTheGroup : groupsOfCellsToFind )
    {
        VectorLong &cellsOfTheGroupToSend = groupsOfCellsToSend[nameOfTheGroup];
        VectorLong &cellsOfTheGroupGathered = groupsOfCellsGathered[nameOfTheGroup];

        AsterMPI::all_gather( cellsOfTheGroupToSend, cellsOfTheGroupGathered );
        cellsOfTheGroupToSend.clear();
    }

    /* ************************************************************************ */
    /* Part of the code dealing with the definition of ConnectionMesh variables */
    /* ************************************************************************ */

    /* Renumbering to have always the same thing */
    for ( int i = 0; i < totalNumberOfNodes; ++i )
    {
        const int globalNum = numNodesGathered[3 * i + 1];
        renumNodes[globalNum] = i;
    }

    pos = 0;
    VectorLong renumNodesLocNew(totalNumberOfNodes, -1);
    for ( auto& it : renumNodes )
    {
        renumNodesLocNew[it.second] = pos;
        pos++;
    }
    renumNodes.clear();


    /* Add numbering */
    _nodesLocalNumbering->allocate( totalNumberOfNodes );
    _nodesGlobalNumbering->allocate( totalNumberOfNodes );
    _nodesOwner->allocate( totalNumberOfNodes );

    for ( int i = 0; i < totalNumberOfNodes; ++i )
    {
        ( *_nodesLocalNumbering )[renumNodesLocNew[i]] = numNodesGathered[3 * i];
        ( *_nodesGlobalNumbering )[renumNodesLocNew[i]] = numNodesGathered[3 * i + 1];
        ( *_nodesOwner )[renumNodesLocNew[i]] = numNodesGathered[3 * i + 2];
        numNodesGloLoc[ ( *_nodesGlobalNumbering )[renumNodesLocNew[i]] ] = renumNodesLocNew[i] + 1;
    }

    /* Add coordinates */
    const auto numberOfConnectionMeshCoordinates=coordinatesGathered.size();
    *_coordinates->getDescriptor() = *mesh->getCoordinates()->getDescriptor();
    *_coordinates->getReference() = *mesh->getCoordinates()->getReference();
    auto values = _coordinates->getValues();
    values->allocate( numberOfConnectionMeshCoordinates );
    values->updateValuePointer();
    for ( int i = 0; i < totalNumberOfNodes; ++i )
    {
        ( *values )[3*renumNodesLocNew[i]]   = coordinatesGathered[3*i];
        ( *values )[3*renumNodesLocNew[i]+1] = coordinatesGathered[3*i+1];
        ( *values )[3*renumNodesLocNew[i]+2] = coordinatesGathered[3*i+2];
    }

    /* Add nodes */
    _nameOfNodes->allocate( totalNumberOfNodes );
    for ( int i = 0; i < totalNumberOfNodes; ++i )
    {
        std::stringstream sstream;
        const int newNum = renumNodesLocNew[i] + 1;
        sstream << std::setfill( '0' ) << std::setw( 7 ) << std::hex << newNum;
        _nameOfNodes->add( newNum, std::string( "N" + sstream.str() ) );
    }

    /* Add group of nodes */
    if( groupsOfNodesToFind.size() > 0 )
    {
        _groupsOfNodes->allocate( groupsOfNodesToFind.size() );

        for ( const auto &nameOfTheGroup : groupsOfNodesToFind ) {
            const auto &toCopy = groupsOfNodesGathered[nameOfTheGroup];
            const auto nbNodes = toCopy.size();
            VectorLong nodesOfGrp( nbNodes );
            for( int i = 0; i < nbNodes; i++)
            {
                auto it = numNodesGloLoc.find( toCopy[i] );
                if( it == numNodesGloLoc.end() )
                    throw std::runtime_error( "Not finding nodes");

                nodesOfGrp[ i ] = it->second;
            }

            _groupsOfNodes->allocateObjectByName( nameOfTheGroup, nbNodes );
            _groupsOfNodes->getObjectFromName( nameOfTheGroup ).setValues( nodesOfGrp );
        }
    }

    /* Add cells */
    AS_ASSERT(totalNumberOfCells > 0);
    _cellsLocalNumbering->allocate( totalNumberOfCells );
    _cellsOwner->allocate( totalNumberOfCells );
    _nameOfCells->allocate( totalNumberOfCells );
    _cellsType->allocate( totalNumberOfCells );
    _connectivity->allocateContiguous( totalNumberOfCells,
                                         connectivitiesGathered.size(), Numbered );

    int offset = 0;
    for ( int i = 1; i <= totalNumberOfCells; ++i )
    {
        std::stringstream sstream;
        sstream << std::setfill( '0' ) << std::setw( 7 ) << std::hex << i;
        _nameOfCells->add( i, std::string( "M" + sstream.str() ) );
        ( *_cellsType )[i - 1] = connectivitiesGathered[offset++];
        ( *_cellsLocalNumbering )[i - 1] = connectivitiesGathered[offset++];
        const auto globCellId = connectivitiesGathered[offset++];
        numCellsGloLoc[ globCellId ] = i;
        ( *_cellsOwner )[i - 1] = connectivitiesGathered[offset++];

        const auto nbNodes = connectivitiesGathered[offset++];
        VectorLong listNodes( nbNodes );
        for( int iNode = 0; iNode < nbNodes; iNode++)
        {
            const auto globNodeId = connectivitiesGathered[offset++];
            auto it = numNodesGloLoc.find( globNodeId );
                if( it == numNodesGloLoc.end() )
                    throw std::runtime_error( "Not finding nodes");

            listNodes[iNode] = it->second;
        }

        _connectivity->allocateObject( nbNodes );
        _connectivity->getObject( i ).setValues( listNodes );
    }
    connectivitiesGathered.clear();

    /* Add group of nodes */
    if( groupsOfCellsToFind.size() > 0 )
    {
        _groupsOfCells->allocate( groupsOfCellsToFind.size() );

        for ( const auto &nameOfTheGroup : groupsOfCellsToFind ) {
            const auto &toCopy = groupsOfCellsGathered[nameOfTheGroup];
            const auto nbCells = toCopy.size();
            VectorLong cellsOfGrp( nbCells );
            for( int i = 0; i < nbCells; i++)
            {
                auto it = numCellsGloLoc.find( toCopy[i] );
                if( it == numCellsGloLoc.end() )
                    throw std::runtime_error( "Not finding cells");

                if( it->second < 0 )
                    throw std::runtime_error( "Wrong cell index" );

                cellsOfGrp[ i ] = it->second;
            }

            _groupsOfCells->allocateObjectByName( nameOfTheGroup, nbCells );
            _groupsOfCells->getObjectFromName( nameOfTheGroup ).setValues( cellsOfGrp );
        }
    }

    /* Update information about the mesh */
    _dimensionInformations->allocate( 6 );
    ( *_dimensionInformations )[0] = totalNumberOfNodes;
    ( *_dimensionInformations )[2] = totalNumberOfCells;
    ( *_dimensionInformations )[5] = mesh->getDimension();
    CALLO_CARGEO( getName() );

    build();
};

VectorString ConnectionMesh::getGroupsOfCells( ) const {
    ASTERINTEGER size = _nameOfGrpCells->size();
    VectorString names;
    for ( int i = 0; i < size; i++ ) {
        names.push_back( trim( _nameOfGrpCells->getStringFromIndex( i + 1 ) ) );
    }
    return names;
}

VectorString ConnectionMesh::getGroupsOfNodes( ) const {
    ASTERINTEGER size = _nameOfGrpNodes->size();
    VectorString names;
    for ( int i = 0; i < size; i++ ) {
        names.push_back( trim( _nameOfGrpNodes->getStringFromIndex( i + 1 ) ) );
    }
    return names;
}


VectorLong ConnectionMesh::getCellsGlobalNumbering( const JeveuxVectorLong& rankOfCells ) const
{
    /* Local variables gathering the basic MPI informations needed */
    AsterMPI AsterMPI;

    /* Rank of the current MPI proc */
    const int rank = getMPIRank();
    /* Total number of processors */
    const int numberOfProcessors = getMPISize();

    rankOfCells->updateValuePointer();
    const int nbCells = _pMesh->getNumberOfCells();
    int nbCellOwned = 0;
    for( int i = 0; i < nbCells; i++ )
    {
        if( ( *rankOfCells )[i] == rank )
            nbCellOwned++;
    }

    /* Sizes obtained after having MPI gathered all sizeInfo variables */
    VectorInt sizePerRank( numberOfProcessors, -1 );
    /* Cumulated sizes obtained in sizePerRank (starting from zero) */
    VectorInt sizeOffset( numberOfProcessors, -1 );

    /* Gather the number of cell owned on all procs */
    AsterMPI::all_gather( nbCellOwned, sizePerRank );

    sizeOffset[0] = 0;
    for ( int i = 1; i < numberOfProcessors; ++i )
    {
        sizeOffset[i] = sizeOffset[i - 1] + sizePerRank[i - 1];
    }

    VectorLong globalNum( nbCells, -1 );

    int globalId = sizeOffset[ rank ];
    for( int i = 0; i < nbCells; i++ )
    {
        if( ( *rankOfCells )[ i ] == rank )
            globalNum[ i ] = globalId++;
    }

    return globalNum;
};

bool ConnectionMesh::hasGroupOfCells( const std::string &name) const {
    if ( _groupsOfCells->size() < 0 && !_groupsOfCells->build() ) {
        return false;
    }
    return _groupsOfCells->existsObject( name );
}

bool ConnectionMesh::hasGroupOfNodes( const std::string &name) const {
    if ( _groupsOfNodes->size() < 0 && !_groupsOfNodes->build() ) {
        return false;
    }
    return _groupsOfNodes->existsObject( name );
}

VectorLong ConnectionMesh::getCells( const std::string name ) const {

    if ( name.empty())
    {
        return irange(long(1), long(getNumberOfCells()));
    }
    else if ( !hasGroupOfCells( name ) ) {
        return VectorLong();
    }

    return _groupsOfCells->getObjectFromName( name ).toVector();
};



#endif /* ASTER_HAVE_MPI */
