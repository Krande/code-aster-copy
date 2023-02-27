/**
 * @file DOFNumbering.cxx
 * @brief Implementation de DOFNumbering
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

#include "Numbering/BaseDOFNumbering.h"

#include "aster_fort_calcul.h"
#include "astercxx.h"

#include "Supervis/ResultNaming.h"

BaseDOFNumbering::BaseDOFNumbering( const std::string name, const std::string &type )
    : DataStructure( name, 14, type ),
      _nameOfSolverDataStructure( JeveuxVectorChar24( getName() + ".NSLV" ) ),
      _smos( new MorseStorage( getName() + ".SMOS" ) ),
      _slcs( new LigneDeCiel( getName() + ".SLCS" ) ),
      _mltf( new MultFrontGarbage( getName() + ".MLTF" ) ),
      _localNumbering( std::make_shared< LocalEquationNumbering >( getName() ) ){};

bool BaseDOFNumbering::computeNumbering( const std::vector< MatrElem > matrix ) {
    ASTERINTEGER nb_matr = matrix.size();
    JeveuxVectorChar24 jvListOfMatr( ResultNaming::getNewResultName() );
    jvListOfMatr->allocate( nb_matr );

    int ind = 0;
    BaseMeshPtr savedMesh( nullptr );
    for ( const auto &mat : matrix ) {
        ( *jvListOfMatr )[ind++] = std::visit( ElementaryMatrixGetName(), mat );
        auto FEDescs = std::visit( ElementaryMatrixGetFEDescrp(), mat );
        this->addFiniteElementDescriptor( FEDescs );
        if ( ind == 1 ) {
            savedMesh = std::visit( ElementaryMatrixGetMesh(), mat );
            if ( savedMesh == nullptr )
                throw std::runtime_error( "No mesh in ElementaryMatrix" );
        } else {
            auto tmp = std::visit( ElementaryMatrixGetMesh(), mat );
            if ( savedMesh != nullptr && tmp != nullptr ) {
                if ( savedMesh->getName() != tmp->getName() )
                    throw std::runtime_error( "Inconsistent mesh in list of ElementaryMatrix" );
            }
        }
    }

    CALLO_NUME_DDL_MATR( getName(), jvListOfMatr->getName(), &nb_matr );

    for ( const auto &mat : matrix ) {
        const auto model = std::visit( ElementaryMatrixGetModel(), mat );
        if ( model ) {
            this->setModel( model );
        } else {
            this->setMesh( std::visit( ElementaryMatrixGetMesh(), mat ) );
        }
    }

    return true;
};

bool BaseDOFNumbering::computeNumbering( const ModelPtr model, const ListOfLoadsPtr listOfLoads ) {
    listOfLoads->build( model );

    const std::string base( "GG" );
    const std::string null( " " );
    CALLO_NUMERO_WRAP( getName(), base, null, null, model->getName(), listOfLoads->getName() );

    const auto FEDescs = listOfLoads->getFiniteElementDescriptors();
    this->addFiniteElementDescriptor( FEDescs );
    this->setModel( model );

    return true;
};

bool BaseDOFNumbering::computeRenumbering( const ModelPtr model,
                                           const ListOfLoadsPtr listOfLoads ) {
    listOfLoads->build( model );

    const std::string base( "GG" );
    const std::string null( " " );

    CALLO_NUMER3( model->getName(), listOfLoads->getName(), getName(), null, base );
    this->setModel( model );

    return true;
};

bool BaseDOFNumbering::computeNumbering( const std::vector< FiniteElementDescriptorPtr > &Feds,
                                         const std::string &localMode ) {

    JeveuxVectorChar24 list_ligrel( "&&LIST_LIGREL" );
    list_ligrel->reserve( Feds.size() );

    for ( auto &fed : Feds ) {
        this->addFiniteElementDescriptor( fed );
        list_ligrel->push_back( JeveuxChar24( fed->getName() ) );
    }

    CALLO_NUME_DDL_CHAMELEM( getName(), list_ligrel->getName(), localMode );

    return true;
};

bool BaseDOFNumbering::addFiniteElementDescriptor( const FiniteElementDescriptorPtr &curFED ) {
    if ( curFED ) {
        const auto name = trim( curFED->getName() );
        if ( _FEDNames.find( name ) == _FEDNames.end() ) {
            _FEDVector.push_back( curFED );
            _FEDNames.insert( name );
            return true;
        }
    }
    return false;
};

bool BaseDOFNumbering::addFiniteElementDescriptor(
    const std::vector< FiniteElementDescriptorPtr > &curFED ) {

    for ( auto &fed : curFED ) {
        auto ret = this->addFiniteElementDescriptor( fed );
        if ( !ret ) {
            return false;
        }
    }

    return true;
};
