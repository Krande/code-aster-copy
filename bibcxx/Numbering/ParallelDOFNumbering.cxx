/**
 * @file ParallelDOFNumbering.cxx
 * @brief Implementation de ParallelDOFNumbering
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

#include "Numbering/ParallelDOFNumbering.h"

#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

#include <stdexcept>

#ifdef ASTER_HAVE_MPI

ParallelDOFNumbering::ParallelDOFNumbering()
    : ParallelDOFNumbering( ResultNaming::getNewResultName() ) {};

ParallelDOFNumbering::ParallelDOFNumbering( const std::string name,
                                            const ParallelEquationNumberingPtr globNume,
                                            const ModelPtr model )
    : BaseDOFNumbering( name, "NUME_DDL_P" ), _globalNumbering( globNume ) {
    setModel( model );
};

ParallelDOFNumbering::ParallelDOFNumbering( const std::string &name )
    : BaseDOFNumbering( name, "NUME_DDL_P" ),
      _globalNumbering( std::make_shared< ParallelEquationNumbering >( getName() + ".NUME" ) ) {};

bool ParallelDOFNumbering::useLagrangeDOF() const { return _globalNumbering->useLagrangeDOF(); };

VectorLong ParallelDOFNumbering::getPhysicalDOFs( const bool local ) const {
    return _globalNumbering->getPhysicalDOFs( local );
};

VectorLong ParallelDOFNumbering::getGhostDOFs( const bool local ) const {
    return _globalNumbering->getGhostDOFs( local );
};

VectorLong ParallelDOFNumbering::getNoGhostDOFs( const bool local ) const {
    return _globalNumbering->getNoGhostDOFs( local );
};

VectorLong ParallelDOFNumbering::getLagrangeDOFs( const bool local ) const {
    return _globalNumbering->getLagrangeDOFs( local );
};

std::map< ASTERINTEGER, VectorLong >
ParallelDOFNumbering::getDictOfLagrangeDOFs( const bool local ) const {
    return _globalNumbering->getDictOfLagrangeDOFs( local );
};

std::string ParallelDOFNumbering::getComponentFromDOF( const ASTERINTEGER dof,
                                                       const bool local ) const {
    return _globalNumbering->getComponentFromDOF( dof, local );
};

ASTERINTEGER
ParallelDOFNumbering::getNodeFromDOF( const ASTERINTEGER dof, const bool local ) const {
    return _globalNumbering->getNodeFromDOF( dof, local );
};

bool ParallelDOFNumbering::isPhysicalDOF( const ASTERINTEGER dof, const bool local ) const {
    return _globalNumbering->isPhysicalDOF( dof, local );
};

ASTERINTEGER
ParallelDOFNumbering::getNumberOfDOFs( const bool local ) const {
    return _globalNumbering->getNumberOfDOFs( local );
};

bool ParallelDOFNumbering::useSingleLagrangeDOF() const {
    return _globalNumbering->useSingleLagrangeDOF();
};

VectorString ParallelDOFNumbering::getComponents() const {
    return _globalNumbering->getComponents();
};

const JeveuxVectorLong ParallelDOFNumbering::getLocalToGlobalMapping() const {
    return _globalNumbering->getLocalToGlobal();
};

ASTERINTEGER ParallelDOFNumbering::localToGlobalDOF( const ASTERINTEGER loc ) const {
    return _globalNumbering->localToGlobalDOF( loc );
}

std::vector< std::pair< ASTERINTEGER, std::string > >
ParallelDOFNumbering::getNodeAndComponentFromDOF( const bool local ) const {
    return _globalNumbering->getNodeAndComponentFromDOF( local );
};

std::pair< ASTERINTEGER, std::string >
ParallelDOFNumbering::getNodeAndComponentFromDOF( const ASTERINTEGER dof, const bool local ) const {
    return _globalNumbering->getNodeAndComponentFromDOF( dof, local );
};

ASTERINTEGER ParallelDOFNumbering::globalToLocalDOF( const ASTERINTEGER glob ) const {
    return _globalNumbering->globalToLocalDOF( glob );
};

VectorString ParallelDOFNumbering::getComponentFromNode( const ASTERINTEGER node,
                                                         const bool local ) const {
    auto localnode = node;
    if ( !local )
        localnode = getMesh()->getGlobalToLocalNodeId( node );
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
        stringVector.push_back( strip( std::string( stringArray + 8 * k, 8 ) ) );
    }
    FreeStr( stringArray );
    return stringVector;
};

bool ParallelDOFNumbering::build() { return _globalNumbering->build(); }

#endif /* ASTER_HAVE_MPI */
