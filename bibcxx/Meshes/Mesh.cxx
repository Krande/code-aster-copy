/**
 * @file Mesh.cxx
 * @brief Implementation de Mesh
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

#include "astercxx.h"

#include "Meshes/Mesh.h"
#include "Utilities/Tools.h"


bool Mesh::readAsterFile( const std::string &fileName ) {
    readMeshFile( fileName, "ASTER" );
    return true;
}

bool Mesh::readGibiFile( const std::string &fileName ) {
    readMeshFile( fileName, "GIBI" );
    return true;
}

bool Mesh::readGmshFile( const std::string &fileName ) {
    readMeshFile( fileName, "GMSH" );
    return true;
}

bool Mesh::hasGroupOfCells( const std::string &name, const bool local ) const {
    if ( _groupsOfCells->size() < 0 && !_groupsOfCells->build() ) {
        return false;
    }
    return _groupsOfCells->existsObject( name );
}

bool Mesh::hasGroupOfNodes( const std::string &name, const bool local ) const {
    if ( _groupsOfNodes->size() < 0 && !_groupsOfNodes->build() ) {
        return false;
    }
    return _groupsOfNodes->existsObject( name );
}

VectorString Mesh::getGroupsOfCells(const bool local) const {
    ASTERINTEGER size = _nameOfGrpCells->size();
    VectorString names;
    for ( int i = 0; i < size; i++ ) {
        names.push_back( trim( _nameOfGrpCells->getStringFromIndex( i + 1 ) ) );
    }
    return names;
}

VectorString Mesh::getGroupsOfNodes(const bool local) const {
    ASTERINTEGER size = _nameOfGrpNodes->size();
    VectorString names;
    for ( int i = 0; i < size; i++ ) {
        names.push_back( trim( _nameOfGrpNodes->getStringFromIndex( i + 1 ) ) );
    }
    return names;
}

VectorLong Mesh::getCells( const std::string name ) const {

    if ( name.empty() ) {
        return irange( long( 1 ), long( getNumberOfCells() ) );
    } else if ( !hasGroupOfCells( name ) ) {
        return VectorLong();
    }

    return _groupsOfCells->getObjectFromName( name ).toVector();
}

VectorLong Mesh::getNodes( const std::string name, const bool localNumbering,
                           const bool same_rank ) const {
    if ( name.empty() ) {
        return irange( long( 1 ), long( getNumberOfNodes() ) );
    } else if ( !hasGroupOfNodes( name ) ) {
        return VectorLong();
    }
    return _groupsOfNodes->getObjectFromName( name ).toVector();
}

VectorLong Mesh::getNodesFromCells( const std::string name, const bool localNumbering,
                                    const bool same_rank ) const {
    CALL_JEMARQ();
    const auto cellsId = getCells( name );

    if(cellsId.empty())
        return VectorLong();

    const auto &connecExp = getConnectivityExplorer();

    SetLong nodes;

    for ( auto &cellId : cellsId ) {
        const auto cell = connecExp[cellId];
        for (auto &node : cell)
            nodes.insert(node);
    }

    CALL_JEDEMA();
    return VectorLong( nodes.begin(), nodes.end() );
};

bool Mesh::isQuadratic() const {
    CALL_JEMARQ();

    auto cellsType = getMedCellsTypes();
    cellsType->updateValuePointer();
    const auto nb_elem = cellsType->size();

    for ( ASTERINTEGER ii = 0; ii < nb_elem; ii++ ) {
        const auto cellType = ( *cellsType )[ii];
        if ( cellType == 103 || cellType == 104 || cellType == 206 || cellType == 207 ||
             cellType == 208 || cellType == 209 || cellType == 310 || cellType == 315 ||
             cellType == 318 || cellType == 313 || cellType == 320 || cellType == 327 )
            return true;
    }

    CALL_JEDEMA();

    return false;
};
