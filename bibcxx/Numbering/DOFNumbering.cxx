/**
 * @file DOFNumbering.cxx
 * @brief Implementation de DOFNumbering
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

#include "astercxx.h"

#include "Numbering/DOFNumbering.h"

#include "Supervis/CommandSyntax.h"
#include "Supervis/ResultNaming.h"

#include <stdexcept>

DOFNumbering::DOFNumbering() : DOFNumbering( ResultNaming::getNewResultName() ) {};

DOFNumbering::DOFNumbering( const std::string name, const EquationNumberingPtr globNume,
                            const ModelPtr model )
    : BaseDOFNumbering( name, "NUME_DDL" ), _globalNumbering( globNume ) {
    setModel( model );
};

DOFNumbering::DOFNumbering( const std::string name )
    : BaseDOFNumbering( name, "NUME_DDL" ),
      _globalNumbering( std::make_shared< EquationNumbering >( getName() + ".NUME" ) ) {};

bool DOFNumbering::useLagrangeDOF() const { return getEquationNumbering()->useLagrangeDOF(); };

bool DOFNumbering::useSingleLagrangeDOF() const {
    return getEquationNumbering()->useSingleLagrangeDOF();
};

VectorLong DOFNumbering::getPhysicalDOFs( const bool local ) const {
    return getEquationNumbering()->getPhysicalDOFs( local );
};

VectorLong DOFNumbering::getNoGhostDOFs( const bool local ) const {
    return getEquationNumbering()->getNoGhostDOFs( local );
};

VectorLong DOFNumbering::getLagrangeDOFs( const bool local ) const {
    return getEquationNumbering()->getLagrangeDOFs( local );
};

std::map< ASTERINTEGER, VectorLong > DOFNumbering::getDictOfLagrangeDOFs( const bool local ) const {
    return getEquationNumbering()->getDictOfLagrangeDOFs( local );
}

std::string DOFNumbering::getComponentFromDOF( const ASTERINTEGER dof, const bool local ) const {
    return getEquationNumbering()->getComponentFromDOF( dof, local );
};

ASTERINTEGER DOFNumbering::getNodeFromDOF( const ASTERINTEGER dof, const bool local ) const {
    return getEquationNumbering()->getNodeFromDOF( dof, local );
};

ASTERINTEGER DOFNumbering::getDOFFromNodeAndComponent( const ASTERINTEGER &node,
                                                       const std::string &comp,
                                                       const bool local ) const {
    return getEquationNumbering()->getDOFFromNodeAndComponent( node, comp, local );
}

std::vector< std::pair< ASTERINTEGER, std::string > >
DOFNumbering::getNodeAndComponentFromDOF( const bool local ) const {
    return _globalNumbering->getNodeAndComponentFromDOF( local );
};

std::pair< ASTERINTEGER, std::string >
DOFNumbering::getNodeAndComponentFromDOF( const ASTERINTEGER dof, const bool local ) const {
    return _globalNumbering->getNodeAndComponentFromDOF( dof, local );
};

bool DOFNumbering::isPhysicalDOF( const ASTERINTEGER dof, const bool local ) const {
    return getEquationNumbering()->isPhysicalDOF( dof, local );
};

ASTERINTEGER DOFNumbering::getNumberOfDOFs( const bool local ) const {
    return getEquationNumbering()->getNumberOfDOFs();
};

VectorString DOFNumbering::getComponents() const {
    return getEquationNumbering()->getComponents();
};

VectorString DOFNumbering::getComponentFromNode( const ASTERINTEGER node, const bool local ) const {
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
        stringVector.push_back( strip( std::string( stringArray + 8 * k, 8 ) ) );
    }
    FreeStr( stringArray );
    return stringVector;
};
