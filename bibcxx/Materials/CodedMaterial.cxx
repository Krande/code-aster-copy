/**
 * @file CodedMaterial.cxx
 * @brief Implementation de CodedMaterial
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

#include "Materials/CodedMaterial.h"
#include "aster_fort_jeveux.h"
#include "aster_fort_material.h"
#include "aster_fort_utils.h"

CodedMaterial::CodedMaterial( const std::string &name, const MaterialFieldPtr &mater,
                                        const ModelPtr &model )
    : DataStructure( name, 19, "MATER_CODE" ), _mater( mater ), _model( model ),
      _field( boost::make_shared< ConstantFieldOnCellsLong >
        ( getName() , _model->getMesh()) ),
      _grp( JeveuxVectorChar8( getName() + ".GRP" ) ),
      _nGrp( JeveuxVectorLong( getName() + ".NGRP" ) ){};

bool CodedMaterial::allocate( bool force ) {
    if ( !force && _field->exists() )
        return false;
    if ( _field->exists() ) {
        _field->deallocate();
        _grp->deallocate();
        _nGrp->deallocate();
        _vecOfCodiVectors.clear();
        _vecOfR8.clear();
        _vecOfIa.clear();
    }
    std::string blanc( 24, ' ' );
    std::string materName = _mater->getName();
    materName.resize( 24, ' ' );
    std::string mate = blanc;
    bool thm( _model->existsThm() );
    // TODO existsTher ?
    bool ther( false );
    std::string strJeveuxBase( "G" );
    CALLO_RCMFMC( materName, mate, (ASTERLOGICAL *)&thm, (ASTERLOGICAL *)&ther, getName(),
                  strJeveuxBase );

    auto vecOfMater = _mater->getVectorOfMaterial();
    for ( auto curIter : vecOfMater ) {
        // Fill codivectors (can be optimized)
        std::string nameWithoutBlanks = trim( getName() ) + ".0";
        std::string base( " " );
        ASTERINTEGER pos = 1;
        ASTERINTEGER nbval2 = 0;
        ASTERINTEGER retour = 0;
        JeveuxChar24 nothing( " " );
        CALLO_JELSTC( base, nameWithoutBlanks, &pos, &nbval2, nothing, &retour );
        if ( retour != 0 ) {
            JeveuxVectorChar24 test( -retour );
            ASTERINTEGER nbval2 = -retour;
            CALLO_JELSTC( base, nameWithoutBlanks, &pos, &nbval2, ( *test )[0], &retour );
            for ( int i = 0; i < retour; ++i ) {
                std::string name = ( *test )[i].toString();
                std::string name2( name, 19, 5 );
                if ( name2 == ".CODI" )
                    _vecOfCodiVectors.push_back( JeveuxVectorLong( name ) );
            }
        }

        const int nbMB = curIter->getNumberOfMaterialProperties();
        for ( int i = 0; i < nbMB; ++i ) {
            auto vecVec1 = curIter->getBehaviourVectorOfValuesReal( i );
            auto vecVec2 = curIter->getBehaviourVectorOfFunctions( i );
            for ( auto vec1 : vecVec1 )
                if ( vec1->exists() ) {
                    const std::string name( vec1->getName(), 0, 16 );
                    const std::string name1 = name + ".LISV_VR";
                    const std::string name2 = name + ".LISV_IA";
                    _vecOfR8.push_back( JeveuxVectorReal( name1 ) );
                    _vecOfIa.push_back( JeveuxVectorLong( name2 ) );
                }
            for ( auto vec2 : vecVec2 )
                if ( vec2->exists() ) {
                    const std::string name( vec2->getName(), 0, 16 );
                    const std::string name1 = name + ".LISV_VR";
                    const std::string name2 = name + ".LISV_IA";
                    _vecOfR8.push_back( JeveuxVectorReal( name1 ) );
                    _vecOfIa.push_back( JeveuxVectorLong( name2 ) );
                }
        }
    }
    return true;
};

bool CodedMaterial::constant() const {
    const std::string typeco( "CHAM_MATER" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "ELAS_FO" );

    CALLO_DISMOI( questi, _mater->getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        return false;
    return true;
};
