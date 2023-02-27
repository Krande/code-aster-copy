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
    : ParallelDOFNumbering( ResultNaming::getNewResultName() ){};

ParallelDOFNumbering::ParallelDOFNumbering( const std::string name,
                                            const ParallelGlobalEquationNumberingPtr globNume,
                                            const ModelPtr model )
    : BaseDOFNumbering( name, "NUME_DDL_P" ), _globalNumbering( globNume ) {
    setModel( model );
};

ParallelDOFNumbering::ParallelDOFNumbering( const std::string &name )
    : BaseDOFNumbering( name, "NUME_DDL_P" ),
      _globalNumbering(
          std::make_shared< ParallelGlobalEquationNumbering >( getName() + ".NUME" ) ){};

bool ParallelDOFNumbering::useLagrangeMultipliers() const {
    return getGlobalEquationNumbering()->useLagrangeMultipliers();
};

VectorLong ParallelDOFNumbering::getRowsAssociatedToPhysicalDofs( const bool local ) const {
    return getGlobalEquationNumbering()->getRowsAssociatedToPhysicalDofs( local );
};

VectorLong ParallelDOFNumbering::getGhostRows( const bool local ) const {
    return getGlobalEquationNumbering()->getGhostRows( local );
};

VectorLong ParallelDOFNumbering::getNoGhostRows() const {
    return getGlobalEquationNumbering()->getNoGhostRows();
};

VectorLong ParallelDOFNumbering::getRowsAssociatedToLagrangeMultipliers( const bool local ) const {
    return getGlobalEquationNumbering()->getRowsAssociatedToLagrangeMultipliers( local );
};

std::string ParallelDOFNumbering::getComponentAssociatedToRow( const ASTERINTEGER row,
                                                               const bool local ) const {
    return getGlobalEquationNumbering()->getComponentAssociatedToRow( row, local );
};

ASTERINTEGER
ParallelDOFNumbering::getNodeAssociatedToRow( const ASTERINTEGER row, const bool local ) const {
    return getGlobalEquationNumbering()->getNodeAssociatedToRow( row, local );
};

bool ParallelDOFNumbering::isRowAssociatedToPhysical( const ASTERINTEGER row,
                                                      const bool local ) const {
    return getGlobalEquationNumbering()->isRowAssociatedToPhysical( row, local );
};

ASTERINTEGER
ParallelDOFNumbering::getNumberOfDofs( const bool local ) const {
    return getGlobalEquationNumbering()->getNumberOfDofs( local );
};

bool ParallelDOFNumbering::useSingleLagrangeMultipliers() const {
    return getGlobalEquationNumbering()->useSingleLagrangeMultipliers();
};

VectorString ParallelDOFNumbering::getComponents() const {
    return getGlobalEquationNumbering()->getComponents();
};

const JeveuxVectorLong ParallelDOFNumbering::getLocalToGlobalMapping() const {
    return getGlobalEquationNumbering()->getLocalToGlobal();
};

const ASTERINTEGER ParallelDOFNumbering::localToGlobalRow( const ASTERINTEGER loc ) {
    return getGlobalEquationNumbering()->localToGlobalRow( loc );
}

const ASTERINTEGER ParallelDOFNumbering::globalToLocalRow( const ASTERINTEGER glob ) const {
    return getGlobalEquationNumbering()->globalToLocalRow( glob );
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

#endif /* ASTER_HAVE_MPI */
