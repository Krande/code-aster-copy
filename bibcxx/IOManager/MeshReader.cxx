/**
 * @file MeshReader.cxx
 * @brief Implementation de MeshReader
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

#include "IOManager/MeshReader.h"

#include "IOManager/MedFileReader.h"

#include <algorithm>

#ifdef ASTER_HAVE_MED

static const VectorInt asterTypeList = { 1,  2,  4,  6,  7,  9,  11, 12, 14, 16,
                                         18, 19, 20, 21, 22, 23, 24, 25, 26, 27 };

static const std::map< int, med_int > asterMedMatching = {
    { 1, 1 },    { 2, 102 },  { 4, 103 },  { 6, 104 },  { 7, 203 },  { 9, 206 },  { 11, 207 },
    { 12, 204 }, { 14, 208 }, { 16, 209 }, { 18, 304 }, { 19, 310 }, { 20, 306 }, { 21, 315 },
    { 22, 318 }, { 23, 305 }, { 24, 313 }, { 25, 308 }, { 26, 320 }, { 27, 327 }
};

const std::set< med_int > medTypeToRenumber = { 304, 308, 305, 306, 310, 320, 313, 315, 318, 327 };

template < std::size_t N, const int ( &indices )[N] >
void applyPermutation( const ASTERINTEGER *in, ASTERINTEGER *out ) {
    for ( size_t i = 0; i < N; i++ ) {
        out[i] = in[indices[i]];
    }
};

VectorLong medToAsterRenumbering( const med_int &medType, const VectorLong &toRenumber,
                                  const int &nbElem ) {
    VectorLong out( toRenumber.size() );
    switch ( medType ) {
    case 304: {
        constexpr std::size_t N = 4;
        static constexpr int arr[N] { 0, 2, 1, 3 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 308: {
        constexpr std::size_t N = 8;
        static constexpr int arr[N] { 0, 3, 2, 1, 4, 7, 6, 5 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 305: {
        constexpr std::size_t N = 5;
        static constexpr int arr[N] { 0, 3, 2, 1, 4 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 306: {
        constexpr std::size_t N = 6;
        static constexpr int arr[N] { 0, 2, 1, 3, 5, 4 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 310: {
        constexpr std::size_t N = 10;
        static constexpr int arr[N] { 0, 2, 1, 3, 6, 5, 4, 7, 9, 8 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 320: {
        constexpr std::size_t N = 20;
        static constexpr int arr[N] { 0, 3, 2,  1,  4,  7,  6,  5,  11, 10,
                                      9, 8, 16, 19, 18, 17, 15, 14, 13, 12 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 313: {
        constexpr std::size_t N = 13;
        static constexpr int arr[N] { 0, 3, 2, 1, 4, 8, 7, 6, 5, 9, 12, 11, 10 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 315: {
        constexpr std::size_t N = 15;
        static constexpr int arr[N] { 0, 2, 1, 3, 5, 4, 8, 7, 6, 12, 14, 13, 11, 10, 9 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 318: {
        constexpr std::size_t N = 18;
        static constexpr int arr[N] {
            0, 2, 1, 3, 5, 4, 8, 7, 6, 12, 14, 13, 11, 10, 9, 17, 16, 15
        };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 327: {
        constexpr std::size_t N = 27;
        static constexpr int arr[N] { 0,  3,  2,  1,  4,  7,  6,  5,  11, 10, 9,  8,  16, 19,
                                      18, 17, 15, 14, 13, 12, 20, 24, 23, 22, 21, 25, 26 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    }
    return out;
}

MeshPtr MeshReader::readFromMedFile( const std::filesystem::path &filename ) {
    MeshPtr toReturn = std::make_shared< Mesh >();
    auto coordsToFill = toReturn->getCoordinates();
    auto coordValues = coordsToFill->getValues();

    // Read mesh from file
    auto fr = MedFileReader();
    fr.openParallel( filename, MedReadOnly );
    const auto curMesh = fr.getMesh( 0 );
    const auto seq = curMesh->getSequence( 0 );
    const auto nodeNb = curMesh->getNodeNumberAtSequence( seq[0], seq[1] );

    // Read node coordinates and copy in coordValues
    const auto nbNodes = nodeNb;
    const auto curCoords = curMesh->readCoordinates( seq[0], seq[1] );
    const auto dim = curMesh->getDimension();
    if ( dim == 3 ) {
        *( coordValues ) = curMesh->readCoordinates( seq[0], seq[1] );
    } else {
        throw std::runtime_error( "Not yet implemented" );
    }

    // Get cell type and sort it according to aster sort
    auto cellTypes = curMesh->getGeometricTypesAtSequence( seq[0], seq[1] );
    std::set< med_int > typeSet( cellTypes.begin(), cellTypes.end() );
    std::vector< med_int > cellTypesSorted;
    VectorInt asterCellTypes;
    for ( const auto &asterType : asterTypeList ) {
        const auto &medType = asterMedMatching.at( asterType );
        if ( typeSet.count( medType ) != 0 ) {
            cellTypesSorted.push_back( medType );
            asterCellTypes.push_back( asterType );
        }
    }

    // Get cell informations by type
    auto connectivity = toReturn->getConnectivity();
    auto cellType = toReturn->getCellTypeVector();
    int totalSize = 0, size = 0;
    VectorInt elemNbAndSizeVec;
    for ( const auto medType : cellTypesSorted ) {
        const auto cellNb =
            curMesh->getCellNumberForGeometricTypeAtSequence( seq[0], seq[1], medType );
        const auto nbNodesForGeoT = curMesh->getNodeNumberForGeometricType( medType );
        totalSize += cellNb * nbNodesForGeoT;
        size += cellNb;
        elemNbAndSizeVec.push_back( cellNb );
        elemNbAndSizeVec.push_back( nbNodesForGeoT );
    }

    // Get families in mesh
    const auto families = curMesh->getFamilies();
    med_int maxId = 0, minId = 0, nodeGrpCount = 0, cellGrpCount = 0;
    std::map< int, MedFamilyPtr > idToFamily;
    VectorString nodeGroupList, cellGroupList;
    std::map< std::string, int > nodeGroupNameToGroupId, cellGroupNameToGroupId;
    for ( const auto &fam : families ) {
        maxId = std::max( fam->getId(), maxId );
        minId = std::min( fam->getId(), minId );
        idToFamily[fam->getId()] = fam;
        const auto &curGroups = fam->getGroups();
        for ( const auto &grpName : curGroups ) {
            if ( fam->getId() > 0 ) {
                if ( nodeGroupNameToGroupId.count( grpName ) == 0 ) {
                    nodeGroupNameToGroupId[grpName] = nodeGrpCount;
                    nodeGroupList.push_back( grpName );
                    ++nodeGrpCount;
                }

            } else {
                if ( cellGroupNameToGroupId.count( grpName ) == 0 ) {
                    cellGroupNameToGroupId[grpName] = cellGrpCount;
                    cellGroupList.push_back( grpName );
                    ++cellGrpCount;
                }
            }
        }
    }
    const int familyOffset = -minId;
    maxId += familyOffset;
    VectorOfVectorsLong cellFamily( maxId + 1, VectorLong() );

    // Get node families
    auto curNFam = curMesh->getNodeFamilyAtSequence( seq[0], seq[1] );
    for ( const auto &[index, cellFamId] : enumerate( curNFam ) ) {
        if ( cellFamId != 0 ) {
            cellFamily[cellFamId + familyOffset].push_back( index + 1 );
        }
    }

    // allocate connectivity
    connectivity->allocate( size, totalSize );
    cellType->allocate( size );
    int count = 0, cumElem = 1, totalCount = 0;
    ASTERINTEGER *connexPtr = nullptr;
    for ( const auto medType : cellTypesSorted ) {
        const auto &cellNb = elemNbAndSizeVec[count * 2];
        const auto &nbNodesForGeoT = elemNbAndSizeVec[count * 2 + 1];

        // Collection allocation
        for ( int cellId = 0; cellId < cellNb; ++cellId ) {
            connectivity->fastAllocateObject( cumElem + cellId, nbNodesForGeoT );
            if ( cellId % 1000000 == 0 )
                std::cout << "Boucle 3 " << medType << " " << cumElem + cellId << std::endl
                          << std::flush;
        }
        // Get contiguous collection start
        if ( count == 0 ) {
            connexPtr = connectivity->startPointer();
        }

        // Read connectivity from file
        auto curConn =
            curMesh->getConnectivityForGeometricTypeAtSequence( seq[0], seq[1], medType );
        // Apply renumbering
        if ( medTypeToRenumber.count( medType ) != 0 ) {
            curConn = medToAsterRenumbering( medType, curConn, cellNb );
        }
        // And copy to collection
        std::copy( curConn.begin(), curConn.end(), connexPtr + totalCount );
        curConn.clear();

        // Fill cell type
        auto cellTypePtr = &( ( *cellType )[cumElem - 1] );
        const auto &curAsterType = asterCellTypes[count];
        std::fill( cellTypePtr, cellTypePtr + cellNb, curAsterType );

        // Get cell families
        auto curFam = curMesh->getCellFamilyForGeometricTypeAtSequence( seq[0], seq[1], medType );
        for ( const auto &[index, cellFamId] : enumerate( curFam ) ) {
            if ( cellFamId != 0 ) {
                cellFamily[cellFamId + familyOffset].push_back( index + cumElem );
            }
        }

        cumElem += cellNb;
        totalCount += cellNb * nbNodesForGeoT;
        ++count;
    }

    int index = 0;
    VectorOfVectorsLong nodeIdGroupList( nodeGroupList.size(), VectorLong() );
    VectorOfVectorsLong cellIdGroupList( cellGroupList.size(), VectorLong() );
    // From families to groups
    for ( const auto &cellFamVector : cellFamily ) {
        const auto famId = index - familyOffset;
        if ( cellFamVector.size() != 0 ) {
            const auto &curFam = idToFamily.at( famId );
            const auto &curGroups = curFam->getGroups();
            for ( const auto &group : curGroups ) {
                if ( famId > 0 ) {
                    const auto &grpId = nodeGroupNameToGroupId.at( group );
                    auto &curGrpToFill = nodeIdGroupList[grpId];
                    for ( const auto &id : cellFamVector ) {
                        curGrpToFill.push_back( id );
                    }
                } else if ( famId < 0 ) {
                    const auto &grpId = cellGroupNameToGroupId.at( group );
                    auto &curGrpToFill = cellIdGroupList[grpId];
                    for ( const auto &id : cellFamVector ) {
                        curGrpToFill.push_back( id );
                    }
                }
            }
        }
        ++index;
    }
    if ( nodeGroupList.size() != 0 ) {
        toReturn->addGroupsOfNodes( nodeGroupList, nodeIdGroupList );
    }
    if ( cellGroupList.size() != 0 ) {
        toReturn->addGroupsOfCells( cellGroupList, cellIdGroupList );
    }
    toReturn->buildInformations( dim );
    toReturn->buildNamesVectors();
    toReturn->endDefinition();
    return toReturn;
}

// void add_automatic_names2( NamesMapChar8 &map, int size, std::string prefix ) {
//     map->allocate( size );
//     if ( size > 10000000 ) {
//         for ( auto i = 0; i < size; ++i ) {
//             std::ostringstream oss;
//             oss << std::hex << i + 1;
//             std::string name = prefix + toUpper( oss.str() );
//             map->add( i + 1, name );
//         }
//     } else {
//         for ( auto i = 0; i < size; ++i ) {
//             map->add( i + 1, prefix + std::to_string( i + 1 ) );
//         }
//     }
// }

// void MeshReader::testPerf() {
//     NamesMapChar8 _nameOfNodes( NamesMapChar8( "TOTO.NOMNOE    " ) );
//     auto nbCells = 55000000;
//     std::cout << "nbCells " << nbCells << std::endl;
//     add_automatic_names2( _nameOfNodes, nbCells, "N" );
// }

#endif
