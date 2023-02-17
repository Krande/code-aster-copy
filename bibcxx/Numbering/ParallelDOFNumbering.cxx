/**
 * @file ParallelDOFNumbering.cxx
 * @brief Implementation de ParallelDOFNumbering
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

#include "Numbering/ParallelDOFNumbering.h"

#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

#include <stdexcept>

#ifdef ASTER_HAVE_MPI

ParallelDOFNumbering::ParallelDOFNumbering()
    : BaseDOFNumbering( ResultNaming::getNewResultName(), "NUME_DDL_P" ) {
    _globalNumbering = std::make_shared< ParallelGlobalEquationNumbering >( getName() );
};

ParallelDOFNumbering::ParallelDOFNumbering( const std::string name,
                                            const FieldOnNodesDescriptionPtr fdof,
                                            const ModelPtr model )
    : BaseDOFNumbering( name, "NUME_DDL_P", fdof ) {
    _globalNumbering = std::make_shared< ParallelGlobalEquationNumbering >( getName() );
    setModel( model );
};

ParallelDOFNumbering::ParallelDOFNumbering( const std::string &name )
    : BaseDOFNumbering( name, "NUME_DDL_P" ) {
    _globalNumbering = std::make_shared< ParallelGlobalEquationNumbering >( getName() );
};

bool ParallelDOFNumbering::useLagrangeMultipliers() const {
    const std::string typeco( "NUME_DDL" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "EXIS_LAGR" );
    bool local_answer = false, global_answer;

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        local_answer = true;

    AsterMPI::all_reduce( local_answer, global_answer, MPI_LOR );

    return global_answer;
};

VectorLong ParallelDOFNumbering::getRowsAssociatedToPhysicalDofs( const bool local ) const {
    auto dofInformation = getGlobalNumbering()->getLagrangianInformations();
    dofInformation->updateValuePointer();
    ASTERINTEGER size = dofInformation->size();
    VectorLong physicalRows;
    ASTERINTEGER physicalIndicator;
    auto loc2glo = getLocalToGlobalMapping();
    loc2glo->updateValuePointer();

    for ( int i = 0; i < size; i++ ) {
        physicalIndicator = ( *dofInformation )[i];
        if ( physicalIndicator == 0 )
            if ( local )
                physicalRows.push_back( i );
            else {
                physicalRows.push_back( ( *loc2glo )[i] );
            }
    }
    return physicalRows;
};

VectorLong ParallelDOFNumbering::getGhostRows( const bool local ) const {
    auto localToRank = getGlobalNumbering()->getLocalToRank();
    localToRank->updateValuePointer();
    VectorLong ghostRows;
    const auto rank = getMPIRank();
    ASTERINTEGER dofOwner;
    auto loc2glo = getLocalToGlobalMapping();
    if ( !local )
        loc2glo->updateValuePointer();

    for ( int i = 0; i < getNumberOfDofs( true ); i++ ) {
        dofOwner = ( *localToRank )[i];
        if ( dofOwner != rank )
            if ( local )
                ghostRows.push_back( i );
            else {
                ghostRows.push_back( ( *loc2glo )[i] );
            }
    }
    return ghostRows;
};

VectorLong ParallelDOFNumbering::getNoGhostRows() const {
    auto localToRank = getGlobalNumbering()->getLocalToRank();
    localToRank->updateValuePointer();
    const auto rank = getMPIRank();
    ASTERINTEGER dofOwner;
    VectorLong noGhostRows;

    for ( int i = 0; i < getNumberOfDofs( true ); i++ ) {
        dofOwner = ( *localToRank )[i];
        if ( dofOwner == rank )
            noGhostRows.push_back( i );
    }
    return noGhostRows;
};

VectorLong ParallelDOFNumbering::getRowsAssociatedToLagrangeMultipliers( const bool local ) const {
    auto dofInformation = getGlobalNumbering()->getLagrangianInformations();
    dofInformation->updateValuePointer();
    ASTERINTEGER size = dofInformation->size();
    VectorLong lagrangeRows;
    ASTERINTEGER physicalIndicator;
    auto loc2glo = getLocalToGlobalMapping();
    loc2glo->updateValuePointer();

    for ( int i = 0; i < size; i++ ) {
        physicalIndicator = ( *dofInformation )[i];
        if ( physicalIndicator != 0 )
            if ( local )
                lagrangeRows.push_back( i );
            else {
                lagrangeRows.push_back( ( *loc2glo )[i] );
            }
    }
    return lagrangeRows;
};

std::string ParallelDOFNumbering::getComponentAssociatedToRow( const ASTERINTEGER row,
                                                               const bool local ) const {
    auto localrow = row;
    if ( !local )
        localrow = globalToLocalRow( row );

    auto [nodeId, cmpName] = getDescription()->getNodeAndComponentFromDOF( localrow );
    return cmpName;
};

ASTERINTEGER
ParallelDOFNumbering::getRowAssociatedToNodeComponent( const ASTERINTEGER node,
                                                       const std::string compoName,
                                                       const bool local ) const {
    auto localnode = node;
    auto loc2glo = getLocalToGlobalMapping();
    loc2glo->updateValuePointer();
    if ( !local )
        localnode =
            std::static_pointer_cast< ParallelMesh >( getMesh() )->globalToLocalNodeId( node );
    if ( localnode < 0 or localnode >= getMesh()->getNumberOfNodes() )
        throw std::out_of_range( "Invalid node index" );
    NamesMapChar8 nodeNameMap = getMesh()->getNameOfNodesMap();
    const std::string nodeName = nodeNameMap->getStringFromIndex( localnode + 1 );
    const std::string objectType( "NUME_DDL" );
    ASTERINTEGER node2, row;

    CALLO_POSDDL( objectType, getName(), nodeName, compoName, &node2, &row );
    assert( localnode + 1 == node2 );
    if ( node2 == 0 )
        throw std::out_of_range( "No node " + std::to_string( node2 ) + " in the mesh" );
    if ( row == 0 )
        throw std::runtime_error( "Node " + std::to_string( node2 ) + " has no " + compoName +
                                  " dof" );
    auto outrow = local ? row - 1 : ( *loc2glo )[row - 1];
    return outrow;
};

ASTERINTEGER
ParallelDOFNumbering::getNodeAssociatedToRow( const ASTERINTEGER row, const bool local ) const {
    auto localrow = row;
    if ( !local )
        localrow = globalToLocalRow( row );

    auto [nodeId, cmpId] = getDescription()->getNodeAndComponentNumberFromDOF( localrow, local );

    return nodeId;
};

bool ParallelDOFNumbering::isRowAssociatedToPhysical( const ASTERINTEGER row,
                                                      const bool local ) const {
    auto localrow = row;
    if ( !local )
        localrow = globalToLocalRow( row );
    auto [nodeId, cmpId] = getDescription()->getNodeAndComponentNumberFromDOF( localrow );

    return cmpId > 0;
};

ASTERINTEGER
ParallelDOFNumbering::getNumberOfDofs( const bool local ) const {
    getGlobalNumbering()->getNumberOfEquations()->updateValuePointer();
    if ( local )
        return ( *getGlobalNumbering()->getNumberOfEquations() )[0];
    else
        return ( *getGlobalNumbering()->getNumberOfEquations() )[1];
};

bool ParallelDOFNumbering::useSingleLagrangeMultipliers() const {
    const std::string typeco( "NUME_DDL" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "SIMP_LAGR" );
    bool local_answer = false, global_answer;

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        local_answer = true;
    AsterMPI::all_reduce( local_answer, global_answer, MPI_LAND );

    return global_answer;
};

VectorString ParallelDOFNumbering::getComponents() const {
    ASTERINTEGER ncmp, maxCmp = 100, ibid = 0;
    char *stringArray;
    VectorString localComp, globalComp;
    std::string all( "ALL" );
    stringArray = MakeTabFStr( 8, maxCmp );
    CALL_NUMEDDL_GET_COMPONENTS( getName().c_str(), all.c_str(), &ibid, &ncmp, stringArray,
                                 &maxCmp );
    for ( int k = 0; k < ncmp; k++ ) {
        localComp.push_back( trim( std::string( stringArray + 8 * k, 8 ) ) );
    }
    FreeStr( stringArray );

    // Communicate with others
    AsterMPI::all_gather( localComp, globalComp );

    return unique( globalComp );
};

VectorString ParallelDOFNumbering::getComponentsAssociatedToNode( const ASTERINTEGER node,
                                                                  const bool local ) const {
    auto localnode = node;
    if ( !local )
        localnode =
            std::static_pointer_cast< ParallelMesh >( getMesh() )->globalToLocalNodeId( node );
    if ( localnode < 0 or localnode >= getMesh()->getNumberOfNodes() )
        throw std::out_of_range( "Invalid node index" );
    ASTERINTEGER ncmp, maxCmp = 100;
    char *stringArray;
    VectorString stringVector;
    std::string all( "ONE" );
    stringArray = MakeTabFStr( 8, maxCmp );
    if ( node < 0 or node >= getMesh()->getNumberOfNodes() )
        throw std::out_of_range( "Invalid node index" );
    ASTERINTEGER aster_node = node + 1;
    CALL_NUMEDDL_GET_COMPONENTS( getName().c_str(), all.c_str(), &aster_node, &ncmp, stringArray,
                                 &maxCmp );
    for ( int k = 0; k < ncmp; k++ ) {
        stringVector.push_back( trim( std::string( stringArray + 8 * k, 8 ) ) );
    }
    FreeStr( stringArray );
    return stringVector;
};

const JeveuxVectorLong ParallelDOFNumbering::getLocalToGlobalMapping() const {
    return getGlobalNumbering()->getLocalToGlobal();
};

const ASTERINTEGER ParallelDOFNumbering::localToGlobalRow( const ASTERINTEGER loc ) {
    if ( loc < 0 or loc >= getNumberOfDofs( true ) )
        throw std::out_of_range( "Invalid row index" );
    auto loc2glo = getLocalToGlobalMapping();
    loc2glo->updateValuePointer();
    return ( *loc2glo )[loc];
}

const ASTERINTEGER ParallelDOFNumbering::globalToLocalRow( const ASTERINTEGER glob ) const {
    if ( _global2localMap.empty() )
        const_cast< ParallelDOFNumbering * >( this )->_buildGlobal2LocalMap();

    auto search = _global2localMap.find( glob );
    if ( search != _global2localMap.end() ) {
        return search->second;
    } else {
        auto rank = getMPIRank();
        throw std::out_of_range( "Global Dof number " + std::to_string( glob ) +
                                 " not found on rank " + std::to_string( rank ) );
        return -1;
    }
};

#endif /* ASTER_HAVE_MPI */
