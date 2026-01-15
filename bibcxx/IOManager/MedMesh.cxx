/**
 * @file MedMesh.cxx
 * @brief Implementation de MedMesh
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
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

// aslint: disable=C3010

#include "IOManager/MedMesh.h"

#include "IOManager/MedFilter.h"
#include "IOManager/MedTypes.h"
#include "IOManager/MedUtilities.h"
#include "ParallelUtilities/AsterMPI.h"

#ifdef ASTER_HAVE_MED
MedMesh::MedMesh( const MedFilePointer &filePtr, const std::string &name, med_int dim )
    : _name( name ), _filePtr( filePtr ), _dim( dim ) {};

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
    auto nbElemT = getCellNumberAtSequence( numdt, numit, iterator );
    auto nbNoCell = getNodeNumberForGeometricType( geotype );
    std::vector< med_int > connectivity;
    if ( _filePtr.isParallel() ) {
        const auto rank = getMPIRank();
        const auto nbProcs = getMPISize();
        const auto pair = splitEntitySet( nbElemT, rank, nbProcs );
        int nbElemL = pair.first;
        int start = pair.second;
        MedFilter medFilter( _filePtr, nbElemT, 1, nbNoCell, MED_ALL_CONSTITUENT,
                             MED_FULL_INTERLACE, MED_COMPACT_STMODE, start, nbElemL, 1, nbElemL,
                             0 );
        connectivity = std::vector< med_int >( nbElemL * nbNoCell, 0 );
        MEDmeshElementConnectivityAdvancedRd( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                                              MED_CELL, geotype, MED_NODAL, medFilter.getPointer(),
                                              &connectivity[0] );
    } else {
        connectivity = std::vector< med_int >( nbElemT * nbNoCell, 0 );
        MEDmeshElementConnectivityRd( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_CELL,
                                      geotype, MED_NODAL, MED_FULL_INTERLACE, &connectivity[0] );
    }
    return connectivity;
};

std::vector< med_int >
MedMesh::getConnectivityForGeometricTypeAtSequence( int numdt, int numit,
                                                    med_geometry_type geotype ) const {
    auto nbElemT = getCellNumberForGeometricTypeAtSequence( numdt, numit, geotype );
    auto nbNoCell = getNodeNumberForGeometricType( geotype );
    std::vector< med_int > connectivity;
    if ( _filePtr.isParallel() ) {
        const auto rank = getMPIRank();
        const auto nbProcs = getMPISize();
        const auto pair = splitEntitySet( nbElemT, rank, nbProcs );
        int nbElemL = pair.first;
        int start = pair.second;
        MedFilter medFilter( _filePtr, nbElemT, 1, nbNoCell, MED_ALL_CONSTITUENT,
                             MED_FULL_INTERLACE, MED_COMPACT_STMODE, start, nbElemL, 1, nbElemL,
                             0 );
        connectivity = std::vector< med_int >( nbElemL * nbNoCell, 0 );
        MEDmeshElementConnectivityAdvancedRd( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                                              MED_CELL, geotype, MED_NODAL, medFilter.getPointer(),
                                              &connectivity[0] );
    } else {
        connectivity = std::vector< med_int >( nbElemT * nbNoCell, 0 );
        MEDmeshElementConnectivityRd( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_CELL,
                                      geotype, MED_NODAL, MED_FULL_INTERLACE, &connectivity[0] );
    }
    return connectivity;
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

std::pair< med_int, med_int >
MedMesh::getSplitCellNumberForGeometricTypeAtSequence( int numdt, int numit,
                                                       med_geometry_type geotype ) const {
    med_bool changement, transformation;
    auto nbMaT =
        MEDmeshnEntity( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_CELL, geotype,
                        MED_CONNECTIVITY, MED_NODAL, &changement, &transformation );
    if ( _filePtr.isParallel() ) {
        const auto rank = getMPIRank();
        const auto nbProcs = getMPISize();
        return splitEntitySet( nbMaT, rank, nbProcs );
    } else {
        return { nbMaT, 1 };
    }
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

std::pair< med_int, med_int > MedMesh::getSplitNodeNumberAtSequence( int numdt, int numit ) const {
    med_bool changement, transformation;
    auto nbNoT =
        MEDmeshnEntity( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_NODE, MED_NO_GEOTYPE,
                        MED_COORDINATE, MED_NO_CMODE, &changement, &transformation );
    if ( _filePtr.isParallel() ) {
        const auto rank = getMPIRank();
        const auto nbProcs = getMPISize();
        return splitEntitySet( nbNoT, rank, nbProcs );
    } else {
        return { nbNoT, 1 };
    }
};

med_int MedMesh::getNodeNumberForGeometricType( med_int geoType ) const {
    med_geometry_type geotype = geoType;
    med_int nbNodes = 0, geoDim = 0;
    MEDmeshGeotypeParameter( _filePtr.getFileId(), geotype, &geoDim, &nbNodes );
    return nbNodes;
}

std::vector< double > MedMesh::readCoordinates( int numdt, int numit ) const {
    std::vector< double > toReturn;
    med_bool changement, transformation;
    auto nbNoT =
        MEDmeshnEntity( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_NODE, MED_NO_GEOTYPE,
                        MED_COORDINATE, MED_NO_CMODE, &changement, &transformation );
    if ( _filePtr.isParallel() ) {
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
        toReturn = std::vector< double >( _dim * nbNoT, 0. );
        MEDmeshNodeCoordinateRd( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                                 MED_FULL_INTERLACE, &toReturn[0] );
    }
    return toReturn;
}

std::vector< med_int > MedMesh::getGlobalNodeNumberingAtSequence( int numdt, int numit ) const {
    med_bool changement, transformation;
    auto nbNoT =
        MEDmeshnEntity( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_NODE, MED_NO_GEOTYPE,
                        MED_COORDINATE, MED_NO_CMODE, &changement, &transformation );
    std::vector< med_int > globNum( nbNoT, 0 );
    MEDmeshGlobalNumberRd( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_NODE,
                           MED_NO_GEOTYPE, &globNum[0] );
    if ( _filePtr.isParallel() ) {
        const auto rank = getMPIRank();
        const auto nbProcs = getMPISize();
        const auto pair = splitEntitySet( nbNoT, rank, nbProcs );
        int nbNoL = pair.first;
        int start = pair.second;
        return std::vector< med_int >( globNum.begin() + start - 1,
                                       globNum.begin() + start - 1 + nbNoL );
    } else {
        return globNum;
    }
};

bool MedMesh::readFromFile() {
    if ( !_filePtr.isOpen() ) {
        return false;
    }

    for ( int j = 1; j <= _nbstep; ++j ) {
        med_int numdt, numit;
        med_float dt;
        MEDmeshComputationStepInfo( _filePtr.getFileId(), _name.c_str(), j, &numdt, &numit, &dt );
        _sequences.emplace_back( numdt, numit, dt );
    }

    const auto famNumber = MEDnFamily( _filePtr.getFileId(), _name.c_str() );
    char *groupname;
    char familyname[MED_NAME_SIZE + 1] = "";
    med_int familynumber;
    for ( int i = 1; i <= famNumber; ++i ) {
        const auto nbGrp = MEDnFamilyGroup( _filePtr.getFileId(), _name.c_str(), i );
        groupname = (char *)malloc( sizeof( char ) * MED_LNAME_SIZE * nbGrp + 1 );
        const auto cret = MEDfamilyInfo( _filePtr.getFileId(), _name.c_str(), i, familyname,
                                         &familynumber, groupname );
        if ( cret < 0 ) {
            const auto natt = MEDnFamily23Attribute( _filePtr.getFileId(), _name.c_str(), i );
            med_int *attval, *attide;
            char *attdes;
            attide = (med_int *)malloc( sizeof( med_int ) * natt );
            attval = (med_int *)malloc( sizeof( med_int ) * natt );
            attdes = (char *)malloc( MED_COMMENT_SIZE * natt + 1 );
            med_int attributenumber = 0, attributevalue = 0;
            MEDfamily23Info( _filePtr.getFileId(), _name.c_str(), i, familyname, attide, attval,
                             attdes, &familynumber, groupname );
            free( attide );
            free( attval );
            free( attdes );
        }
        const auto gnames = splitChar( groupname, nbGrp, MED_LNAME_SIZE );
        _families.emplace_back( new MedFamily( std::string( familyname, strlen( familyname ) ),
                                               familynumber, gnames ) );
        free( groupname );
    }

    const auto jointNumber = MEDnSubdomainJoint( _filePtr.getFileId(), _name.c_str() );
    for ( int jointId = 1; jointId <= jointNumber; ++jointId ) {
        char jointname[MED_NAME_SIZE + 1] = "", remotemeshname[MED_NAME_SIZE + 1] = "";
        char description[MED_COMMENT_SIZE + 1] = "";
        med_int domainnumber = -1, nstep = -1, nocstpncorrespondence = -1;
        MEDsubdomainJointInfo( _filePtr.getFileId(), _name.c_str(), jointId, jointname, description,
                               &domainnumber, remotemeshname, &nstep, &nocstpncorrespondence );

        _joints.emplace_back( new MedJoint( _filePtr, _name, jointId, jointname, description,
                                            domainnumber, remotemeshname, nstep,
                                            nocstpncorrespondence ) );
    }
    return true;
};

void MedMesh::printCoordinatesAtSequence( med_int numdt, med_int numit, med_int nodeNumber,
                                          const std::vector< med_float > &coord,
                                          MedFilterPtr filter ) const {
    if ( filter == nullptr ) {
        MEDmeshNodeCoordinateWr( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_NO_DT,
                                 MED_FULL_INTERLACE, nodeNumber, &coord[0] );
    } else {
        MEDmeshNodeCoordinateAdvancedWr( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                                         MED_NO_DT, filter->getPointer(), &coord[0] );
    }
};

void MedMesh::addFamily( const std::string &name, med_int num, const VectorString &grps ) {
    if ( name.size() > MED_NAME_SIZE ) {
        throw std::runtime_error(
            "Family name too long (size limite: " + std::to_string( MED_NAME_SIZE ) + ")" );
    }
    _families.emplace_back( new MedFamily( name, num, grps ) );
    char *grpChar = stringVectorToChar( grps, MED_LNAME_SIZE );
    MEDfamilyCr( _filePtr.getFileId(), _name.c_str(), name.c_str(), num, grps.size(), grpChar );
    free( grpChar );
};

void MedMesh::setNodeFamilyAtSequence( med_int numdt, med_int numit,
                                       const std::vector< med_int > &family, MedFilterPtr filter ) {
    if ( filter == nullptr ) {
        MEDmeshEntityFamilyNumberWr( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_NODE,
                                     MED_NO_GEOTYPE, family.size(), &family[0] );
    } else {
        MEDmeshEntityAttributeAdvancedWr( _filePtr.getFileId(), _name.c_str(), MED_FAMILY_NUMBER,
                                          numdt, numit, MED_NODE, MED_NO_GEOTYPE,
                                          filter->getPointer(), &family[0] );
    }
};

void MedMesh::setCellFamilyAtSequence( med_int numdt, med_int numit, med_geometry_type geotype,
                                       const std::vector< med_int > &family, MedFilterPtr filter ) {
    if ( filter == nullptr ) {
        MEDmeshEntityFamilyNumberWr( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_CELL,
                                     geotype, family.size(), &family[0] );
    } else {
        MEDmeshEntityAttributeAdvancedWr( _filePtr.getFileId(), _name.c_str(), MED_FAMILY_NUMBER,
                                          numdt, numit, MED_CELL, geotype, filter->getPointer(),
                                          &family[0] );
    }
};

void MedMesh::printConnectivityForGeometricTypeAtSequence( med_int numdt, med_int numit,
                                                           med_geometry_type geotype,
                                                           med_int cellNumber,
                                                           const std::vector< med_int > &conn,
                                                           MedFilterPtr filter ) {
    if ( filter == nullptr ) {
        MEDmeshElementConnectivityWr( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_NO_DT,
                                      MED_CELL, geotype, MED_NODAL, MED_FULL_INTERLACE, cellNumber,
                                      &conn[0] );
    } else {
        MEDmeshElementConnectivityAdvancedWr( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                                              MED_NO_DT, MED_CELL, geotype, MED_NODAL,
                                              filter->getPointer(), &conn[0] );
    }
};

void MedMesh::printGlobalNodeNumberingAtSequence( med_int numdt, med_int numit,
                                                  const std::vector< med_int > &globnum ) {
    MEDmeshGlobalNumberWr( _filePtr.getFileId(), _name.c_str(), numdt, numit, MED_NODE,
                           MED_NO_GEOTYPE, globnum.size(), &globnum[0] );
};

void MedMesh::createJoint( const std::string &name, const std::string &desc,
                           med_int oppositeDomain ) {
    if ( name.size() > MED_NAME_SIZE ) {
        throw std::runtime_error(
            "Joint name too long (size limite: " + std::to_string( MED_NAME_SIZE ) + ")" );
    }
    if ( desc.size() > MED_COMMENT_SIZE ) {
        throw std::runtime_error(
            "FDescription too long (size limite: " + std::to_string( MED_COMMENT_SIZE ) + ")" );
    }
    MEDsubdomainJointCr( _filePtr.getFileId(), _name.c_str(), name.c_str(), desc.c_str(),
                         oppositeDomain, _name.c_str() );
};

void MedMesh::printNodeJointAtSequence( const std::string &name, med_int numdt, med_int numit,
                                        const std::vector< med_int > &joint ) {
    MEDsubdomainCorrespondenceWr( _filePtr.getFileId(), _name.c_str(), name.c_str(), numdt, numit,
                                  MED_NODE, MED_NO_GEOTYPE, MED_NODE, MED_NO_GEOTYPE,
                                  joint.size() / 2, &joint[0] );
};

#endif
