/**
 * @file MeshCoordinatesField.cxx
 * @brief Implementation de MeshCoordinatesField
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

#include "DataFields/MeshCoordinatesField.h"

#include "DataFields/FieldOnNodes.h"

/* person_in_charge: nicolas.sellenet at edf.fr */

MeshCoordinatesField &MeshCoordinatesField::operator+=( const FieldOnNodesReal &rhs ) {

    std::string base( "G" ), cumul( "CUMU" );
    ASTERDOUBLE alpha = 1.;
    MeshCoordinatesField oldCoord = this->duplicate();

    CALLO_VTGPLD( cumul, &alpha, oldCoord.getName(), rhs.getName(), base, getName() );

    updateValuePointers();

    return *this;
}

MeshCoordinatesField &MeshCoordinatesField::operator-=( const FieldOnNodesReal &rhs ) {

    std::string base( "G" ), cumul( "CUMU" );
    ASTERDOUBLE alpha = -1.;
    MeshCoordinatesField oldCoord = this->duplicate();

    CALLO_VTGPLD( cumul, &alpha, oldCoord.getName(), rhs.getName(), base, getName() );

    updateValuePointers();

    return *this;
}

VectorLong MeshCoordinatesField::_getDOFsToUse( const VectorString &list_cmp ) const {
    const bool all_nodes = true;

    const bool all_cmp = list_cmp.empty();
    std::map< ASTERINTEGER, std::string > map_cmp;
    if ( !all_cmp ) {
        for ( auto &name : list_cmp ) {
            auto name_trim = trim( name );
            if ( name_trim == "X" ) {
                map_cmp[0] = name_trim;
            } else if ( name_trim == "Y" ) {
                map_cmp[1] = name_trim;
            } else if ( name_trim == "Z" ) {
                map_cmp[2] = name_trim;
            }
        }
    }

    auto taille = this->size();
    VectorLong dofUsed;
    dofUsed.reserve( taille );

    for ( auto dof = 0; dof < taille; ++dof ) {
        bool l_keep_cmp = all_cmp || map_cmp.count( dof % 3 ) > 0;
        bool l_keep_node = all_nodes;
        if ( l_keep_node && l_keep_cmp )
            dofUsed.push_back( dof );
    }

    return dofUsed;
};
