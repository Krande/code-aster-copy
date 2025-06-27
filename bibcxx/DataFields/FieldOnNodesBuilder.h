/**
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

#pragma once

#include "astercxx.h"

#include "aster_fort_ds.h"

#include "DataFields/FieldConverter.h"
#include "DataFields/FieldOnNodes.h"
#include "DataFields/SimpleFieldOnNodes.h"
#include "Modeling/Model.h"
#include "Utilities/Tools.h"

/**
 * There is a separate file to construct some object to avoid circular
 * inclusions. Theses new methods are binding as contructor.
 */

/**
 * @brief Constructor for empty FieldOnNodes.
 */
template < typename ValueType >
std::shared_ptr< FieldOnNodes< ValueType > > FieldOnNodesPtrBuilder( const BaseMeshPtr mesh,
                                                                     const std::string &quantity,
                                                                     const VectorString &cmps ) {

    SimpleFieldOnNodes< ValueType > cns( mesh, quantity, cmps );

    cns.setValues( 0. );

    return toFieldOnNodes( cns );
};

template < typename ValueType >
std::shared_ptr< FieldOnNodes< ValueType > >
FieldOnNodesPtrBuilder( const BaseMeshPtr mesh, const std::string &quantity,
                        const std::map< std::string, ValueType > &values,
                        const VectorString &groupsOfNodes = {},
                        const VectorString &groupsOfCells = {} ) {

    VectorString cmps;
    for ( auto &[cmp, val] : values ) {
        cmps.push_back( cmp );
    }

    SimpleFieldOnNodes< ValueType > cns( mesh, quantity, cmps );

    VectorLong nodes;

    if ( groupsOfCells.empty() || !groupsOfNodes.empty() ) {
        nodes = mesh->getNodes( groupsOfNodes );
    }

    if ( !groupsOfCells.empty() ) {
        auto nodes2 = mesh->getNodesFromCells( groupsOfCells );
        nodes = set_union( nodes, nodes2 );
    }

    cns.setValues( values, nodes );

    return toFieldOnNodes( cns );
};
