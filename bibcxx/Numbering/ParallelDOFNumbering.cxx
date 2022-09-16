/**
 * @file ParallelDOFNumbering.cxx
 * @brief Implementation de ParallelDOFNumbering
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

#include "Numbering/ParallelDOFNumbering.h"

#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

#include <stdexcept>

#ifdef ASTER_HAVE_MPI

ParallelDOFNumbering::ParallelDOFNumbering()
    : BaseDOFNumbering( ResultNaming::getNewResultName(), "NUME_DDL_P" ),
      _globalNumbering( new ParallelGlobalEquationNumbering( getName() ) ) {};

ParallelDOFNumbering::ParallelDOFNumbering( const std::string& name )
    : BaseDOFNumbering( name, "NUME_DDL_P" ),
      _globalNumbering( new ParallelGlobalEquationNumbering( getName() ) ) {};

std::string ParallelDOFNumbering::getPhysicalQuantity() const {
    _globalNumbering->_informations->updateValuePointer();
    JeveuxChar24 physicalQuantity = ( *_globalNumbering->_informations )[1];
    return physicalQuantity.rstrip();
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

    AsterMPI::all_reduce( local_answer, global_answer, MPI_LAND );

    return global_answer;
};

VectorLong ParallelDOFNumbering::getRowsAssociatedToPhysicalDofs( const bool local ) const {
    if ( !local )
        throw std::runtime_error( "Only local operation is supported" );
    getGlobalNumbering()->getLagrangianInformations()->updateValuePointer();
    ASTERINTEGER size = getGlobalNumbering()->getLagrangianInformations()->size();
    VectorLong physicalRows;
    ASTERINTEGER physicalIndicator;
    for ( int i = 0; i < size; i++ ) {
        physicalIndicator = ( *getGlobalNumbering()->getLagrangianInformations() )[i];
        if ( physicalIndicator == 0 )
            physicalRows.push_back( i );
    }
    return physicalRows;
};

VectorLong ParallelDOFNumbering::getRowsAssociatedToLagrangeMultipliers( const bool local ) const {
    if ( !local )
        throw std::runtime_error( "Only local operation is supported" );
    getGlobalNumbering()->getLagrangianInformations()->updateValuePointer();
    ASTERINTEGER size = getGlobalNumbering()->getLagrangianInformations()->size();
    VectorLong lagrangeRows;
    ASTERINTEGER physicalIndicator;
    for ( int i = 0; i < size; i++ ) {
        physicalIndicator = ( *getGlobalNumbering()->getLagrangianInformations() )[i];
        if ( physicalIndicator != 0 )
            lagrangeRows.push_back( i );
    }
    return lagrangeRows;
};

std::string ParallelDOFNumbering::getComponentAssociatedToRow( const ASTERINTEGER row,
                                                               const bool local ) const {
    if ( !local )
        throw std::runtime_error( "Only local operation is supported" );
    if ( row < 0 or row >= getNumberOfDofs() )
        throw std::runtime_error( "Invalid row index" );
    JeveuxVectorLong descriptor = getDescription()->getNodeAndComponentsNumberFromDOF();
    descriptor->updateValuePointer();
    const ASTERINTEGER cmpId = abs( ( *descriptor )[2 * row + 1] );
    if ( cmpId == 0 )
        return " "; // Lagrange multiplier of a MPC - no component
    JeveuxChar8 cmpName( " " );
    CALLO_NUMEDDL_GET_COMPONENT_NAME( getName(), &cmpId, cmpName );

    return cmpName.rstrip();
};

ASTERINTEGER ParallelDOFNumbering::getRowAssociatedToNodeComponent( const ASTERINTEGER node,
                                                                    const std::string compoName,
                                                                    const bool local ) const {
    if ( !local )
        throw std::runtime_error( "Only local operation is supported" );
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

ASTERINTEGER ParallelDOFNumbering::getNodeAssociatedToRow( const ASTERINTEGER row,
                                                           const bool local ) const {
    if ( !local )
        throw std::runtime_error( "Only local operation is supported" );
    if ( row < 0 or row >= getNumberOfDofs( local ) )
        throw std::runtime_error( "Invalid row index" );
    JeveuxVectorLong descriptor = getDescription()->getNodeAndComponentsNumberFromDOF();
    descriptor->updateValuePointer();
    return ( *descriptor )[2 * row] - 1;
};

bool ParallelDOFNumbering::isRowAssociatedToPhysical( const ASTERINTEGER row,
                                                      const bool local ) const {
    if ( row < 0 or row >= getNumberOfDofs( local ) )
        throw std::runtime_error( "Invalid row index" );
    JeveuxVectorLong descriptor = getDescription()->getNodeAndComponentsNumberFromDOF();
    descriptor->updateValuePointer();
    return ( *descriptor )[2 * row + 1] > 0;
};

ASTERINTEGER ParallelDOFNumbering::getNumberOfDofs( const bool local ) const {
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
    if ( !local )
        throw std::runtime_error( "Only local operation is supported" );
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

void ParallelDOFNumbering::setElementaryMatrix(
    const ElementaryMatrixDisplacementRealPtr &currentMatrix ) {
    if ( !currentMatrix->getModel()->getMesh()->isParallel() )
        throw std::runtime_error( "Mesh must be parallel" );
    BaseDOFNumbering::setElementaryMatrix( currentMatrix );
};

void ParallelDOFNumbering::setElementaryMatrix(
    const ElementaryMatrixDisplacementComplexPtr &currentMatrix ) {
    if ( !currentMatrix->getModel()->getMesh()->isParallel() )
        throw std::runtime_error( "Mesh must be parallel" );
    BaseDOFNumbering::setElementaryMatrix( currentMatrix );
};

void ParallelDOFNumbering::setElementaryMatrix(
    const ElementaryMatrixTemperatureRealPtr &currentMatrix ) {
    if ( !currentMatrix->getModel()->getMesh()->isParallel() )
        throw std::runtime_error( "Mesh must be parallel" );
    BaseDOFNumbering::setElementaryMatrix( currentMatrix );
};

void ParallelDOFNumbering::setElementaryMatrix(
    const ElementaryMatrixPressureComplexPtr &currentMatrix ) {
    if ( !currentMatrix->getModel()->getMesh()->isParallel() )
        throw std::runtime_error( "Mesh must be parallel" );
    BaseDOFNumbering::setElementaryMatrix( currentMatrix );
};

void ParallelDOFNumbering::setModel( const ModelPtr &currentModel ) {
    if ( !currentModel->getMesh()->isParallel() )
        throw std::runtime_error( "Mesh must be parallel" );
    BaseDOFNumbering::setModel( currentModel );
};

#endif /* ASTER_HAVE_MPI */
