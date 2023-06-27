/**
 * @file MedMesh.cxx
 * @brief Implementation de MedMesh
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

// aslint: disable=C3010

#include "IOManager/MedMesh.h"

#include "IOManager/MedFilter.h"
#include "IOManager/MedTypes.h"
#include "IOManager/MedUtilities.h"
#include "ParallelUtilities/AsterMPI.h"

MedMesh::MedMesh( const MedFilePointer &filePtr, const std::string &name, med_int dim )
    : _name( name ), _filePtr( filePtr ), _dim( dim ) {
    const auto famNumber = MEDnFamily( _filePtr.getFileId(), _name.c_str() );
    char *groupname;
    for ( int i = 1; i <= famNumber; ++i ) {
        const auto nbGrp = MEDnFamilyGroup( _filePtr.getFileId(), _name.c_str(), i );
        groupname = (char *)malloc( ( nbGrp * MED_LNAME_SIZE + 1 ) * sizeof( char ) );
        char familyname[MED_NAME_SIZE + 1] = "";
        med_int familynumber;
        MEDfamilyInfo( _filePtr.getFileId(), _name.c_str(), i, familyname, &familynumber,
                       groupname );
        const auto gnames = splitChar( groupname, nbGrp, MED_LNAME_SIZE );
        _families.emplace_back( new MedFamily( std::string( familyname, strlen( familyname ) ),
                                               familynumber, gnames ) );
        free( groupname );
    }
};

std::vector< med_int > MedMesh::getCellFamilyAtSequence( int numdt, int numit,
                                                         int iterator ) const {
    med_geometry_type geotype = getCellTypeAtSequence( numdt, numit, iterator );
    auto nbElemT = getCellNumberAtSequence( numdt, numit, iterator );
    std::vector< med_int > family( nbElemT, 0 );
    auto out = MEDmeshEntityFamilyNumberRd( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                                            MED_CELL, geotype, &family[0] );
    if ( _filePtr.isParallel() ) {
        const auto rank = getMPIRank();
        const auto nbProcs = getMPISize();
        const auto pair = splitEntitySet( nbElemT, rank, nbProcs );
        int nbElemL = pair.first;
        int start = pair.second;
        return std::vector< med_int >( family.begin() + start - 1,
                                       family.begin() + start - 1 + nbElemL );
    } else {
        return family;
    }
};

std::vector< med_int >
MedMesh::getCellFamilyForGeometricTypeAtSequence( int numdt, int numit,
                                                  med_geometry_type geotype ) const {
    auto nbElemT = getCellNumberForGeometricTypeAtSequence( numdt, numit, geotype );
    std::vector< med_int > family( nbElemT, 0 );
    auto out = MEDmeshEntityFamilyNumberRd( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                                            MED_CELL, geotype, &family[0] );
    if ( _filePtr.isParallel() ) {
        const auto rank = getMPIRank();
        const auto nbProcs = getMPISize();
        const auto pair = splitEntitySet( nbElemT, rank, nbProcs );
        int nbElemL = pair.first;
        int start = pair.second;
        return std::vector< med_int >( family.begin() + start - 1,
                                       family.begin() + start - 1 + nbElemL );
    } else {
        return family;
    }
};

std::vector< med_int > MedMesh::getNodeFamilyAtSequence( int numdt, int numit ) const {
    med_bool changement, transformation;
    auto nbNoT =
        MEDmeshnEntity( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_NODE, MED_NO_GEOTYPE,
                        MED_COORDINATE, MED_NO_CMODE, &changement, &transformation );
    std::vector< med_int > family( nbNoT, 0 );
    MEDmeshEntityFamilyNumberRd( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_NODE,
                                 MED_NO_GEOTYPE, &family[0] );
    if ( _filePtr.isParallel() ) {
        const auto rank = getMPIRank();
        const auto nbProcs = getMPISize();
        const auto pair = splitEntitySet( nbNoT, rank, nbProcs );
        int nbNoL = pair.first;
        int start = pair.second;
        return std::vector< med_int >( family.begin() + start - 1,
                                       family.begin() + start - 1 + nbNoL );
    } else {
        return family;
    }
};

std::vector< med_int > MedMesh::getConnectivityAtSequence( int numdt, int numit,
                                                           int iterator ) const {
    char geotypename[MED_NAME_SIZE + 1] = "";
    med_geometry_type geotype = MED_NONE;
    MEDmeshEntityInfo( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_CELL, iterator,
                       geotypename, &geotype );
    if ( _filePtr.isParallel() ) {
        auto nbElemT = getCellNumberAtSequence( numdt, numit, iterator );
        auto nbNoCell = getNodeNumberForGeometricType( geotype );
        const auto rank = getMPIRank();
        const auto nbProcs = getMPISize();
        const auto pair = splitEntitySet( nbElemT, rank, nbProcs );
        int nbElemL = pair.first;
        int start = pair.second;
        MedFilter medFilter( _filePtr, nbElemT, 1, nbNoCell, MED_ALL_CONSTITUENT,
                             MED_FULL_INTERLACE, MED_COMPACT_STMODE, start, nbElemL, 1, nbElemL,
                             0 );
        std::vector< med_int > connectivity( nbElemL * nbNoCell, 0 );
        MEDmeshElementConnectivityAdvancedRd( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                                              MED_CELL, geotype, MED_NODAL, medFilter.getPointer(),
                                              &connectivity[0] );
        return connectivity;
    } else {
        throw std::runtime_error( "Not yet implemented" );
    }
};

std::vector< med_int >
MedMesh::getConnectivityForGeometricTypeAtSequence( int numdt, int numit,
                                                    med_geometry_type geotype ) const {
    if ( _filePtr.isParallel() ) {
        auto nbElemT = getCellNumberForGeometricTypeAtSequence( numdt, numit, geotype );
        auto nbNoCell = getNodeNumberForGeometricType( geotype );
        const auto rank = getMPIRank();
        const auto nbProcs = getMPISize();
        const auto pair = splitEntitySet( nbElemT, rank, nbProcs );
        int nbElemL = pair.first;
        int start = pair.second;
        MedFilter medFilter( _filePtr, nbElemT, 1, nbNoCell, MED_ALL_CONSTITUENT,
                             MED_FULL_INTERLACE, MED_COMPACT_STMODE, start, nbElemL, 1, nbElemL,
                             0 );
        std::vector< med_int > connectivity( nbElemL * nbNoCell, 0 );
        MEDmeshElementConnectivityAdvancedRd( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                                              MED_CELL, geotype, MED_NODAL, medFilter.getPointer(),
                                              &connectivity[0] );
        return connectivity;
    } else {
        throw std::runtime_error( "Not yet implemented" );
    }
};

med_int MedMesh::getCellNumberAtSequence( int numdt, int numit, int iterator ) const {
    med_geometry_type geotype = getCellTypeAtSequence( numdt, numit, iterator );
    return getCellNumberForGeometricTypeAtSequence( numdt, numit, geotype );
};

med_int MedMesh::getAllCellNumberAtSequence( int numdt, int numit ) const {
    med_int nbCells = 0;
    for ( int i = 0; i < medTypes.size(); ++i ) {
        const med_geometry_type geotype = medTypes[i];
        const auto test = getCellNumberForGeometricTypeAtSequence( numdt, numit, geotype );
        nbCells += test;
    }
    return nbCells;
};

med_int MedMesh::getCellNumberForGeometricTypeAtSequence( int numdt, int numit,
                                                          med_geometry_type geotype ) const {
    med_bool changement, transformation;
    return MEDmeshnEntity( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_CELL, geotype,
                           MED_CONNECTIVITY, MED_NODAL, &changement, &transformation );
};

med_geometry_type MedMesh::getCellTypeAtSequence( int numdt, int numit, int iterator ) const {
    char geotypename[MED_NAME_SIZE + 1] = "";
    med_geometry_type geotype = MED_NONE;
    MEDmeshEntityInfo( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_CELL, iterator,
                       geotypename, &geotype );
    return geotype;
};

med_int MedMesh::getCellTypeNumberAtSequence( int numdt, int numit ) const {
    med_bool changement, transformation;
    return MEDmeshnEntity( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_CELL, MED_GEO_ALL,
                           MED_UNDEF_DATATYPE, MED_NO_CMODE, &changement, &transformation );
};

std::vector< med_int > MedMesh::getGeometricTypesAtSequence( int numdt, int numit ) const {
    const auto nbCellType = getCellTypeNumberAtSequence( numdt, numit );
    std::vector< med_int > toReturn;
    for ( int i = 0; i < nbCellType; ++i ) {
        toReturn.push_back( getCellTypeAtSequence( numdt, numit, i + 1 ) );
    }
    return toReturn;
}

med_int MedMesh::getNodeNumberAtSequence( int numdt, int numit ) const {
    med_bool changement, transformation;
    return MEDmeshnEntity( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_NODE,
                           MED_NO_GEOTYPE, MED_COORDINATE, MED_NO_CMODE, &changement,
                           &transformation );
};

med_int MedMesh::getNodeNumberForGeometricType( med_int geoType ) const {
    med_geometry_type geotype = geoType;
    med_int nbNodes = 0, geoDim = 0;
    MEDmeshGeotypeParameter( _filePtr.getFileId(), geotype, &geoDim, &nbNodes );
    return nbNodes;
}

std::vector< double > MedMesh::readCoordinates( int numdt, int numit ) const {
    std::vector< double > toReturn;
    if ( _filePtr.isParallel() ) {
        med_bool changement, transformation;
        auto nbNoT = MEDmeshnEntity( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_NODE,
                                     MED_NO_GEOTYPE, MED_COORDINATE, MED_NO_CMODE, &changement,
                                     &transformation );
        const auto rank = getMPIRank();
        const auto nbProcs = getMPISize();
        const auto pair = splitEntitySet( nbNoT, rank, nbProcs );
        int nbNoL = pair.first;
        int start = pair.second;
        MedFilter medFilter( _filePtr, nbNoT, 1, _dim, MED_ALL_CONSTITUENT, MED_FULL_INTERLACE,
                             MED_COMPACT_STMODE, start, nbNoL, 1, nbNoL, 0 );
        toReturn = std::vector< double >( _dim * nbNoL, 0. );
        MEDmeshNodeCoordinateAdvancedRd( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                                         medFilter.getPointer(), &toReturn[0] );
    } else {
        throw std::runtime_error( "Not yet implemented" );
    }
    return toReturn;
}
