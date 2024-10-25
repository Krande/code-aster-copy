/**
 * @file MedProfle.cxx
 * @brief Implementation de MedField
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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
// aslint: disable=C3008

#include "IOManager/MedField.h"

#include "IOManager/MedFilter.h"
#include "IOManager/MedUtilities.h"
#include "ParallelUtilities/AsterMPI.h"

#ifdef ASTER_HAVE_MED
std::vector< med_int > MedField::getAllSupportEntitiesAtSequence( int numdt, int numit ) const {
    std::vector< med_int > toReturn;
    char defaultprofilename[MED_NAME_SIZE + 1] = "";
    char defaultlocalizationname[MED_NAME_SIZE + 1] = "";
    const med_entity_type entitypeN = MED_NODE;
    const med_geometry_type geotypeN = MED_NONE;
    auto nbValN = MEDfieldnProfile( _filePtr.getFileId(), _name.c_str(), numdt, numit, entitypeN,
                                    geotypeN, defaultprofilename, defaultlocalizationname );
    if ( nbValN != 0 ) {
        toReturn.push_back( entitypeN );
        toReturn.push_back( geotypeN );
    }
    const med_entity_type entitype = MED_CELL;
    for ( int i = 0; i < medTypes.size(); ++i ) {
        char defaultprofilename2[MED_NAME_SIZE + 1] = "";
        char defaultlocalizationname2[MED_NAME_SIZE + 1] = "";
        const med_geometry_type geotype = medTypes[i];
        auto nbVal = MEDfieldnProfile( _filePtr.getFileId(), _name.c_str(), numdt, numit, entitype,
                                       geotype, defaultprofilename2, defaultlocalizationname2 );
        if ( nbVal != 0 ) {
            toReturn.push_back( entitype );
            toReturn.push_back( geotype );
        }
    }
    const med_entity_type entitypeNC = MED_NODE_ELEMENT;
    for ( int i = 0; i < medTypes.size(); ++i ) {
        char defaultprofilename2[MED_NAME_SIZE + 1] = "";
        char defaultlocalizationname2[MED_NAME_SIZE + 1] = "";
        const med_geometry_type geotype = medTypes[i];
        auto nbVal =
            MEDfieldnProfile( _filePtr.getFileId(), _name.c_str(), numdt, numit, entitypeNC,
                              geotype, defaultprofilename2, defaultlocalizationname2 );
        if ( nbVal != 0 ) {
            toReturn.push_back( entitypeNC );
            toReturn.push_back( geotype );
        }
    }
    return toReturn;
}

int MedField::getProfileNumberAtSequenceOnEntity( int numdt, int numit, int entitype,
                                                  med_geometry_type geotype ) const {
    char defaultprofilename[MED_NAME_SIZE + 1] = "",
                                            defaultlocalizationname[MED_NAME_SIZE + 1] = "";
    return MEDfieldnProfile( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                             (med_entity_type)entitype, geotype, defaultprofilename,
                             defaultlocalizationname );
};

MedVectorPtr MedField::getValuesAtSequenceOnCellTypesList(
    int numdt, int numit, std::vector< med_geometry_type > medTypesList ) const {
    const bool isParallel = _filePtr.isParallel();

    const med_entity_type entitype = MED_CELL;
    const med_entity_type entitypes[2] = { MED_CELL, MED_NODE_ELEMENT };
    const int csit = _numdtNumitToSeq.at( numdt ).at( numit );
    med_int numdt2, numit2, meshnumdt, meshnumit;
    med_float dt;
    MEDfieldComputingStepMeshInfo( _filePtr.getFileId(), _name.c_str(), csit, &numdt2, &numit2, &dt,
                                   &meshnumdt, &meshnumit );
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();
    int cumulatedElem = 0;
    std::vector< int > nbElems, starts;
    for ( int i = 0; i < medTypesList.size(); ++i ) {
        const med_geometry_type geotype = medTypesList[i];
        const auto nbElemT =
            _mesh->getCellNumberForGeometricTypeAtSequence( meshnumdt, meshnumit, geotype );
        if ( isParallel ) {
            const auto pair = splitEntitySet( nbElemT, rank, nbProcs );
            cumulatedElem += pair.first;
            nbElems.push_back( pair.first );
            starts.push_back( pair.second );
        } else {
            cumulatedElem += nbElemT;
            nbElems.push_back( nbElemT );
            starts.push_back( 1 );
        }
    }
    MedVectorPtr medV( new MedVector( cumulatedElem ) );
    medV->setComponentNumber( _nbCmp );
    medV->setComponentName( _componentName );

    cumulatedElem = 0;
    std::vector< std::vector< med_int > > profs;
    std::vector< med_geometry_type > geoTypes;
    std::vector< std::pair< int, int > > indirection;
    std::vector< int > totSizes;
    std::vector< int > compNb;
    std::vector< int > cumPos;
    std::vector< MedProfilePtr > profsPtr;
    std::vector< int > nbElems2, starts2;
    int verif = 0, medType = -1;
    for ( int entIndex = 0; entIndex < 2; ++entIndex ) {
        const med_entity_type entitype = entitypes[entIndex];
        bool foundSomething = false;
        cumulatedElem = 0;
        for ( int i = 0; i < medTypesList.size(); ++i ) {
            char profilename[MED_NAME_SIZE + 1] = "";
            char localizationname[MED_NAME_SIZE + 1] = "";
            const med_geometry_type geotype = medTypesList[i];
            auto nbProf = MEDfieldnProfile( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                                            entitype, geotype, profilename, localizationname );

            const auto nbCells = nbElems[i];

            if ( nbCells != 0 ) {
                med_int profilesize = 0, nintegrationpoint = 0;
                for ( int profileit = 1; profileit <= nbProf; ++profileit ) {
                    const auto nbVal = MEDfieldnValueWithProfileByName(
                        _filePtr.getFileId(), _name.c_str(), numdt, numit, entitype, geotype,
                        profilename, MED_GLOBAL_STMODE, &profilesize, localizationname,
                        &nintegrationpoint );
                    foundSomething = true;
                    std::string profName( profilename );
                    geoTypes.push_back( geotype );
                    compNb.push_back( _nbCmp );
                    compNb.push_back( nintegrationpoint );
                    totSizes.push_back( nbVal );
                    cumPos.push_back( cumulatedElem );
                    nbElems2.push_back( nbCells );
                    starts2.push_back( starts[i] );
                    if ( profName == "" ) {
                        if ( nbProf != 1 )
                            throw std::runtime_error( "Unexpected value" );
                        medV->setElements( cumulatedElem, nbCells, _nbCmp * nintegrationpoint );
                        profs.emplace_back( 0 );
                        indirection.push_back( { nbElems[i], starts[i] } );
                        profsPtr.push_back( MedProfilePtr( nullptr ) );
                    } else {
                        const auto profSize =
                            MEDprofileSizeByName( _filePtr.getFileId(), profilename );
                        auto &profilearray = profs.emplace_back( profSize );
                        MEDprofileRd( _filePtr.getFileId(), profilename, &profilearray[0] );
                        int deb = -1, taille = 0;
                        const auto start = starts[i];
                        const auto end = nbElems[i] + start;
                        for ( int j = 0; j < profSize; ++j ) {
                            const int curElem = profilearray[j];
                            if ( curElem >= start && curElem < end ) {
                                if ( deb == -1 )
                                    deb = j + 1;
                                medV->setElement( curElem - start + cumulatedElem,
                                                  _nbCmp * nintegrationpoint );
                                ++taille;
                            }
                        }
                        indirection.push_back( { taille, deb } );
                        profsPtr.push_back( _profiles[_mapProfileNameRank.at( profName )] );
                    }
                }
            }
            cumulatedElem += nbCells;
        }
        if ( foundSomething ) {
            medType = entIndex;
            ++verif;
        }
    }
    // There must only be MED_CELL or MED_NODE_ELEMENT but not the two
    // A field is ELGA or ELNO but not the two.
    if ( verif != 1 )
        throw std::runtime_error( "Error while reading med field on cells " + _name );
    medV->endDefinition();
    for ( int i = 0; i < profs.size(); ++i ) {
        const auto &nbNoT = totSizes[i];
        const auto &nbCmp = compNb[2 * i];
        const auto &nbPg = compNb[2 * i + 1];
        const auto &start2 = indirection[i].second;
        const auto &start = starts2[i];
        const auto end = start + nbElems2[i];
        const auto &nbNoL = indirection[i].first;
        const auto &profPtr = profsPtr[i];
        const auto &curPos = cumPos[i];
        const auto &geotype = geoTypes[i];
        const auto cmpPg = nbCmp * nbPg;

        double *value = nullptr;
        VectorReal valVec;
        if ( profPtr == nullptr ) {
            value = &( *medV )[curPos][0];
        } else {
            valVec = VectorReal( cmpPg * nbNoL, 0. );
            value = &valVec[0];
        }
        if ( isParallel ) {
            MedFilter medFilter( _filePtr, nbNoT, nbPg, nbCmp, MED_ALL_CONSTITUENT,
                                 MED_FULL_INTERLACE, MED_COMPACT_STMODE, start2, nbNoL, 1, nbNoL, 0,
                                 profPtr );
            MEDfieldValueAdvancedRd( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                                     entitypes[medType], geotype, medFilter.getPointer(),
                                     (unsigned char *)value );
        } else {
            std::string profilename( "" );
            if ( profPtr != nullptr )
                profilename = profPtr->getName();
            MEDfieldValueWithProfileRd( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                                        entitypes[medType], geotype, MED_COMPACT_STMODE,
                                        profilename.c_str(), MED_FULL_INTERLACE,
                                        MED_ALL_CONSTITUENT, (unsigned char *)value );
        }
        if ( profPtr != nullptr ) {
            const auto &curProf = profs[i];
            const auto &size = curProf.size();
            int count = 0;
            for ( int j = 0; j < size; ++j ) {
                const auto &index = curProf[j];
                if ( index >= start && index < end ) {
                    auto &curElem = ( *medV )[index - start + curPos];
                    for ( int k = 0; k < cmpPg; ++k ) {
                        curElem[k] = value[count];
                        ++count;
                    }
                }
            }
        }
    }
    return medV;
};

MedVectorPtr MedField::getValuesAtSequenceOnNodes( int numdt, int numit ) const {
    const bool isParallel = _filePtr.isParallel();

    const med_entity_type entitype = MED_NODE;
    const int csit = _numdtNumitToSeq.at( numdt ).at( numit );
    med_int numdt2, numit2, meshnumdt, meshnumit;
    med_float dt;
    MEDfieldComputingStepMeshInfo( _filePtr.getFileId(), _name.c_str(), csit, &numdt2, &numit2, &dt,
                                   &meshnumdt, &meshnumit );
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();
    const auto nbNodeT = _mesh->getNodeNumberAtSequence( meshnumdt, meshnumit );
    const std::pair< int, int > pair = ( isParallel ? splitEntitySet( nbNodeT, rank, nbProcs )
                                                    : std::make_pair( (int)nbNodeT, 1 ) );
    const auto &start = pair.second;
    const auto &end = pair.first + start;
    MedVectorPtr medV( new MedVector( pair.first ) );
    medV->setComponentNumber( _nbCmp );
    medV->setComponentName( _componentName );
    medV->setSize( pair.first );

    std::vector< std::vector< med_int > > profs;
    std::vector< std::pair< int, int > > indirection;
    std::vector< int > totSizes;
    std::vector< int > compNb;
    std::vector< MedProfilePtr > profsPtr;

    char profilename[MED_NAME_SIZE + 1] = "";
    char localizationname[MED_NAME_SIZE + 1] = "";
    auto nbProf = MEDfieldnProfile( _filePtr.getFileId(), _name.c_str(), numdt, numit, entitype,
                                    MED_NO_GEOTYPE, profilename, localizationname );

    const auto nbNodes = pair.first;

    if ( nbNodes != 0 ) {
        med_int profilesize = 0, nintegrationpoint = 0;
        for ( int profileit = 1; profileit <= nbProf; ++profileit ) {
            const auto nbVal = MEDfieldnValueWithProfileByName(
                _filePtr.getFileId(), _name.c_str(), numdt, numit, entitype, MED_NO_GEOTYPE,
                profilename, MED_GLOBAL_STMODE, &profilesize, localizationname,
                &nintegrationpoint );
            std::string profName( profilename );
            compNb.push_back( _nbCmp );
            compNb.push_back( nintegrationpoint );
            totSizes.push_back( nbVal );
            if ( profName == "" ) {
                if ( nbProf != 1 )
                    throw std::runtime_error( "Unexpected value" );
                medV->setElements( 0, nbNodes, _nbCmp * nintegrationpoint );
                profs.emplace_back( 0 );
                indirection.push_back( { pair.first, pair.second } );
                profsPtr.push_back( MedProfilePtr( nullptr ) );
            } else {
                const auto profSize = MEDprofileSizeByName( _filePtr.getFileId(), profilename );
                auto &profilearray = profs.emplace_back( profSize );
                MEDprofileRd( _filePtr.getFileId(), profilename, &profilearray[0] );
                int deb = -1, taille = 0;
                for ( int j = 0; j < profSize; ++j ) {
                    const int curElem = profilearray[j];
                    if ( curElem >= start && curElem < end ) {
                        if ( deb == -1 )
                            deb = j + 1;
                        medV->setElement( curElem - start, _nbCmp * nintegrationpoint );
                        ++taille;
                    }
                }
                indirection.push_back( { taille, deb } );
                profsPtr.push_back( _profiles[_mapProfileNameRank.at( profName )] );
            }
        }
    }

    medV->endDefinition();
    for ( int i = 0; i < profs.size(); ++i ) {
        const auto &nbNoT = totSizes[i];
        const auto &nbCmp = compNb[2 * i];
        const auto &nbPg = compNb[2 * i + 1];
        const auto &start2 = indirection[i].second;
        const auto &nbNoL = indirection[i].first;
        const auto &profPtr = profsPtr[i];

        double *value = nullptr;
        VectorReal valVec;
        if ( profPtr == nullptr ) {
            value = &( *medV )[0][0];
        } else {
            valVec = VectorReal( nbCmp * nbPg * nbNoL, 0. );
            value = &valVec[0];
        }
        if ( isParallel ) {
            MedFilter medFilter( _filePtr, nbNoT, nbPg, nbCmp, MED_ALL_CONSTITUENT,
                                 MED_FULL_INTERLACE, MED_COMPACT_STMODE, start2, nbNoL, 1, nbNoL, 0,
                                 profPtr );
            MEDfieldValueAdvancedRd( _filePtr.getFileId(), _name.c_str(), numdt, numit, entitype,
                                     MED_NO_GEOTYPE, medFilter.getPointer(),
                                     (unsigned char *)value );
        } else {
            std::string profilename( "" );
            if ( profPtr != nullptr )
                profilename = profPtr->getName();
            MEDfieldValueWithProfileRd( _filePtr.getFileId(), _name.c_str(), numdt, numit, entitype,
                                        MED_NO_GEOTYPE, MED_COMPACT_STMODE, profilename.c_str(),
                                        MED_FULL_INTERLACE, MED_ALL_CONSTITUENT,
                                        (unsigned char *)value );
        }
        if ( profPtr != nullptr ) {
            const auto &curProf = profs[i];
            const auto &size = curProf.size();
            int count = 0;
            for ( int j = 0; j < size; ++j ) {
                const auto &index = curProf[j];
                if ( index >= start && index < end ) {
                    auto &curElem = ( *medV )[index - start];
                    for ( int k = 0; k < nbCmp * nbPg; ++k ) {
                        curElem[k] = value[count];
                        ++count;
                    }
                }
            }
        }
    }
    return medV;
};

std::vector< double > MedField::getValuesAtSequenceOnEntityAndProfile( int numdt, int numit,
                                                                       int entitype,
                                                                       med_geometry_type geotype,
                                                                       int profileit ) const {
    if ( !_filePtr.isParallel() ) {
        throw std::runtime_error( "Sequential read not yet implemented" );
    }
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();
    char profilename[MED_NAME_SIZE + 1] = "", localizationname[MED_NAME_SIZE + 1] = "";
    med_int profilesize, nintegrationpoint;
    auto nbVal = MEDfieldnValueWithProfile( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                                            (med_entity_type)entitype, geotype, profileit,
                                            MED_GLOBAL_STMODE, profilename, &profilesize,
                                            localizationname, &nintegrationpoint );
    const auto &nbNoT = nbVal;
    const auto &nbCmp = _nbCmp;
    const auto &nbPg = nintegrationpoint;
    const auto &pair = splitEntitySet( nbNoT, rank, nbProcs );
    const auto &start = pair.second;
    const auto &nbNoL = pair.first;
    std::string profName( profilename );
    MedProfilePtr profPtr( nullptr );
    if ( profName != "" ) {
        const auto &profPtr = _profiles[_mapProfileNameRank.at( profName )];
    }
    MedFilter medFilter( _filePtr, nbNoT, nbPg, nbCmp, MED_ALL_CONSTITUENT, MED_FULL_INTERLACE,
                         MED_COMPACT_STMODE, start, nbNoL, 1, nbNoL, 0, profPtr );
    std::vector< double > toReturn( nbNoL * nbCmp * nbPg, 0. );
    MEDfieldValueAdvancedRd( _filePtr.getFileId(), _name.c_str(), numdt, numit,
                             (med_entity_type)entitype, MED_NO_GEOTYPE, medFilter.getPointer(),
                             (unsigned char *)&toReturn[0] );
    return toReturn;
};

#endif
