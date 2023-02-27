/**
 * @file DOFNumbering.cxx
 * @brief Implementation de DOFNumbering
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

#include "Numbering/DOFNumbering.h"

#include "astercxx.h"

#include "Supervis/CommandSyntax.h"
#include "Supervis/ResultNaming.h"

#include <stdexcept>

DOFNumbering::DOFNumbering() : DOFNumbering( ResultNaming::getNewResultName() ){};

DOFNumbering::DOFNumbering( const std::string name, const GlobalEquationNumberingPtr globNume,
                            const ModelPtr model )
    : BaseDOFNumbering( name, "NUME_DDL" ), _globalNumbering( globNume ) {
    setModel( model );
};

DOFNumbering::DOFNumbering( const std::string name )
    : BaseDOFNumbering( name, "NUME_DDL" ),
      _globalNumbering( std::make_shared< GlobalEquationNumbering >( getName() + ".NUME" ) ){};

bool DOFNumbering::useLagrangeMultipliers() const {
    return getGlobalEquationNumbering()->useLagrangeMultipliers();
};

bool DOFNumbering::useSingleLagrangeMultipliers() const {
    return getGlobalEquationNumbering()->useSingleLagrangeMultipliers();
};

VectorLong DOFNumbering::getRowsAssociatedToPhysicalDofs( const bool local ) const {
    return getGlobalEquationNumbering()->getRowsAssociatedToPhysicalDofs( local );
};

VectorLong DOFNumbering::getRowsAssociatedToLagrangeMultipliers( const bool local ) const {
    return getGlobalEquationNumbering()->getRowsAssociatedToLagrangeMultipliers( local );
};

std::string DOFNumbering::getComponentAssociatedToRow( const ASTERINTEGER row,
                                                       const bool local ) const {
    return getGlobalEquationNumbering()->getComponentAssociatedToRow( row, local );
};

ASTERINTEGER DOFNumbering::getNodeAssociatedToRow( const ASTERINTEGER row,
                                                   const bool local ) const {
    return getGlobalEquationNumbering()->getNodeAssociatedToRow( row, local );
};

bool DOFNumbering::isRowAssociatedToPhysical( const ASTERINTEGER row, const bool local ) const {
    return getGlobalEquationNumbering()->isRowAssociatedToPhysical( row, local );
};

ASTERINTEGER DOFNumbering::getNumberOfDofs( const bool local ) const {
    return getGlobalEquationNumbering()->getNumberOfDofs();
};

VectorString DOFNumbering::getComponents() const {
    return getGlobalEquationNumbering()->getComponents();
};

VectorString DOFNumbering::getComponentsAssociatedToNode( const ASTERINTEGER node,
                                                          const bool local ) const {
    ASTERINTEGER ncmp, maxCmp = 100;
    char *stringArray;
    VectorString stringVector;
    std::string all( "ONE" );
    stringArray = MakeTabFStr( 8, maxCmp );
    if ( node < 0 or node >= getMesh()->getNumberOfNodes() )
        throw std::runtime_error( "Invalid node index" );
    ASTERINTEGER aster_node = node + 1;
    CALL_NUMEDDL_GET_COMPONENTS( getName().c_str(), all.c_str(), &aster_node, &ncmp, stringArray,
                                 &maxCmp );
    for ( int k = 0; k < ncmp; k++ ) {
        stringVector.push_back( trim( std::string( stringArray + 8 * k, 8 ) ) );
    }
    FreeStr( stringArray );
    return stringVector;
};

ASTERINTEGER DOFNumbering::getRowAssociatedToNodeComponent( const ASTERINTEGER node,
                                                            const std::string compoName,
                                                            const bool local ) const {
    if ( node < 0 or node >= getMesh()->getNumberOfNodes() )
        throw std::runtime_error( "Invalid node index" );
    NamesMapChar8 nodeNameMap = getMesh()->getNameOfNodesMap();
    const std::string nodeName = nodeNameMap->getStringFromIndex( node + 1 );
    const std::string objectType( "NUME_DDL" );
    ASTERINTEGER node2, row;
    CALLO_POSDDL( objectType, getName(), nodeName, compoName, &node2, &row );
    assert( node + 1 == node2 );
    if ( node2 == 0 )
        throw std::runtime_error( "No node " + std::to_string( node2 ) + " in the mesh" );
    if ( row == 0 )
        throw std::runtime_error( "Node " + std::to_string( node2 ) + " has no " + compoName +
                                  " dof" );
    return row - 1;
};
