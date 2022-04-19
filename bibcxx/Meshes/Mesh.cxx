/**
 * @file Mesh.cxx
 * @brief Implementation de Mesh
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

#include "Meshes/Mesh.h"

#include "astercxx.h"

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

bool Mesh::hasGroupOfCells( const std::string &name, const bool ) const {
    if ( _groupsOfCells->size() < 0 && !_groupsOfCells->build() ) {
        return false;
    }
    return _groupsOfCells->existsObject( name );
}

bool Mesh::hasGroupOfNodes( const std::string &name, const bool ) const {
    if ( _groupsOfNodes->size() < 0 && !_groupsOfNodes->build() ) {
        return false;
    }
    return _groupsOfNodes->existsObject( name );
}

VectorString Mesh::getGroupsOfCells( const bool ) const {
    ASTERINTEGER size = _nameOfGrpCells->size();
    VectorString names;
    for ( auto i = 0; i < size; i++ ) {
        names.push_back( trim( _nameOfGrpCells->getStringFromIndex( i + 1 ) ) );
    }
    return names;
}

VectorString Mesh::getGroupsOfNodes( const bool ) const {
    ASTERINTEGER size = _nameOfGrpNodes->size();
    VectorString names;
    for ( auto i = 0; i < size; i++ ) {
        names.push_back( trim( _nameOfGrpNodes->getStringFromIndex( i + 1 ) ) );
    }
    return names;
}

VectorLong Mesh::getCells( const std::string name ) const {

    if ( name.empty() ) {
        return irange( long( 0 ), long( getNumberOfCells() - 1 ) );
    } else if ( !hasGroupOfCells( name ) ) {
        return VectorLong();
    }
    VectorLong cells = _groupsOfCells->getObjectFromName( name ).toVector();
    for ( auto &cell : cells )
        cell -= 1;
    return cells;
}

VectorLong Mesh::getNodes( const std::string name, const bool, const ASTERINTEGER ) const {
    if ( name.empty() ) {
        return irange( long( 0 ), long( getNumberOfNodes() - 1 ) );
    } else if ( !hasGroupOfNodes( name ) ) {
        return VectorLong();
    }
    VectorLong nodes = _groupsOfNodes->getObjectFromName( name ).toVector();
    for ( auto &node : nodes )
        node -= 1;
    return nodes;
}

VectorLong Mesh::getNodesFromCells( const std::string name, const bool, const ASTERINTEGER ) const {
    CALL_JEMARQ();
    const auto cellsId = getCells( name );

    if ( cellsId.empty() ) {
        CALL_JEDEMA();

        return VectorLong();
    }

    const auto &connecExp = getConnectivityExplorer();

    SetLong nodes;

    for ( auto &cellId : cellsId ) {
        const auto cell = connecExp[cellId];
        for ( auto &node : cell )
            auto ret = nodes.insert( node-1 );
    }

    CALL_JEDEMA();
    return VectorLong( nodes.begin(), nodes.end() );
};

bool Mesh::isQuadratic() const {
    CALL_JEMARQ();

    auto cellsType = getMedCellsTypes();
    cellsType->updateValuePointer();

    for ( auto &cellType : cellsType ) {
        if ( cellType == 103 || cellType == 104 || cellType == 206 || cellType == 207 ||
             cellType == 208 || cellType == 209 || cellType == 310 || cellType == 315 ||
             cellType == 318 || cellType == 313 || cellType == 320 || cellType == 327 ) {
            CALL_JEDEMA();
            return true;
        }
    }

    CALL_JEDEMA();

    return false;
};
