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

DOFNumbering::DOFNumbering()
    : BaseDOFNumbering( ResultNaming::getNewResultName(), "NUME_DDL" ),
      _globalNumbering( std::make_shared< GlobalEquationNumbering >( getName() ) ){};

DOFNumbering::DOFNumbering( const std::string name, const FieldOnNodesDescriptionPtr fdof,
                            const ModelPtr model )
    : BaseDOFNumbering( name, "NUME_DDL", fdof ),
      _globalNumbering( std::make_shared< GlobalEquationNumbering >( getName() ) ) {
    setModel( model );
};

DOFNumbering::DOFNumbering( const std::string name )
    : BaseDOFNumbering( name, "NUME_DDL" ),
      _globalNumbering( std::make_shared< GlobalEquationNumbering >( getName() ) ){};

bool DOFNumbering::useLagrangeMultipliers() const {
    const std::string typeco( "NUME_DDL" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "EXIS_LAGR" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};

VectorLong DOFNumbering::getRowsAssociatedToPhysicalDofs( const bool local ) const {
    auto lagrInfo = getGlobalNumbering()->getLagrangianInformations();
    lagrInfo->updateValuePointer();
    ASTERINTEGER size = lagrInfo->size();
    VectorLong physicalRows;
    ASTERINTEGER physicalIndicator;
    for ( int i = 0; i < size; i++ ) {
        physicalIndicator = ( *lagrInfo )[i];
        if ( physicalIndicator == 0 )
            physicalRows.push_back( i );
    }
    return physicalRows;
};

VectorLong DOFNumbering::getRowsAssociatedToLagrangeMultipliers( const bool local ) const {
    auto lagrInfo = getGlobalNumbering()->getLagrangianInformations();
    lagrInfo->updateValuePointer();
    ASTERINTEGER size = lagrInfo->size();
    VectorLong lagrangeRows;
    ASTERINTEGER physicalIndicator;
    for ( int i = 0; i < size; i++ ) {
        physicalIndicator = ( *lagrInfo )[i];
        if ( physicalIndicator != 0 )
            lagrangeRows.push_back( i );
    }
    return lagrangeRows;
};

std::string DOFNumbering::getComponentAssociatedToRow( const ASTERINTEGER row,
                                                       const bool local ) const {
    auto [nodeId, cmpName] = getDescription()->getNodeAndComponentFromDOF( row );
    return cmpName;
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

ASTERINTEGER DOFNumbering::getNodeAssociatedToRow( const ASTERINTEGER row,
                                                   const bool local ) const {
    auto [nodeId, cmpId] = getDescription()->getNodeAndComponentNumberFromDOF( row );
    return nodeId;
};

bool DOFNumbering::isRowAssociatedToPhysical( const ASTERINTEGER row, const bool local ) const {
    auto [nodeId, cmpId] = getDescription()->getNodeAndComponentNumberFromDOF( row );
    return cmpId > 0;
};

ASTERINTEGER DOFNumbering::getNumberOfDofs( const bool local ) const {
    getGlobalNumbering()->getNumberOfEquations()->updateValuePointer();
    return ( *getGlobalNumbering()->getNumberOfEquations() )[0];
};

bool DOFNumbering::useSingleLagrangeMultipliers() const {
    const std::string typeco( "NUME_DDL" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "SIMP_LAGR" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};

VectorString DOFNumbering::getComponents() const {
    ASTERINTEGER ncmp, maxCmp = 100, ibid = 0;
    char *stringArray;
    VectorString stringVector;
    std::string all( "ALL" );
    stringArray = MakeTabFStr( 8, maxCmp );
    CALL_NUMEDDL_GET_COMPONENTS( getName().c_str(), all.c_str(), &ibid, &ncmp, stringArray,
                                 &maxCmp );
    for ( int k = 0; k < ncmp; k++ ) {
        stringVector.push_back( trim( std::string( stringArray + 8 * k, 8 ) ) );
    }
    FreeStr( stringArray );
    return stringVector;
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
