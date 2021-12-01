/**
 * @file MaterialProperty.cxx
 * @brief Implementation de GenericMaterialProperty
 * @author Nicolas Sellenet
 * @todo autoriser le type Function pour les paramètres matériau
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

#include <stdexcept>
#include "Materials/BaseMaterialProperty.h"

bool GenericMaterialProperty::buildJeveuxVectors(
    JeveuxVectorComplex &complexValues, JeveuxVectorReal &doubleValues,
    JeveuxVectorChar16 &char16Values, JeveuxVectorChar16 &ordr, JeveuxVectorLong &kOrd,
    std::vector< JeveuxVectorReal > &userReals,
    std::vector< JeveuxVectorChar8 > &userFunctions )
{
    const int nbOfMaterialProperties = getNumberOfPropertiesWithValue();
    complexValues->allocate( nbOfMaterialProperties );
    doubleValues->allocate( nbOfMaterialProperties );
    char16Values->allocate( 2 * nbOfMaterialProperties );

    typedef std::map< std::string, int > MapStrInt;
    MapStrInt mapTmp;
    bool bOrdr = false;
    if ( _vectOrdr.size() != 0 ) {
        ordr->allocate( _vectOrdr.size() );
        for ( int i = 0; i < int(_vectOrdr.size()); ++i ) {
            ( *ordr )[i] = _vectOrdr[i];
            mapTmp[_vectOrdr[i]] = i + 1;
        }
        kOrd->allocate( 2 + 4 * nbOfMaterialProperties );
        ( *kOrd )[0] = _vectOrdr.size();
        ( *kOrd )[1] = nbOfMaterialProperties;
        bOrdr = true;
    }

    int position = 0, position2 = nbOfMaterialProperties, pos3 = 0;
    for ( int i = 0; i < int(_vectKW.size()) / 2; ++i ) {
        auto nameOfProperty2 = _vectKW[2 * i];
        auto nameOfProperty = _vectKW[2 * i + 1];
        auto curIter1 = _mapOfRealMaterialProperties.find( nameOfProperty2 );
        if ( curIter1 == _mapOfRealMaterialProperties.end() )
            continue;

        auto curIter = *curIter1;
        if ( curIter.second.hasValue() ) {
            if ( bOrdr ) {
                const int curPos = 2 + nbOfMaterialProperties + i;
                ( *kOrd )[curPos] = position + 1;
                const int curPos2 = 2 + i;
                ( *kOrd )[curPos2] = mapTmp[nameOfProperty];
            }
            nameOfProperty.resize( 16, ' ' );
            ( *char16Values )[position] = nameOfProperty.c_str();
            ( *doubleValues )[position] = curIter.second.getValue();
            ++position;
        }

        if ( curIter.second.isMandatory() && !curIter.second.hasValue() )
            throw std::runtime_error( "Mandatory material property " + nameOfProperty2 +
                                      " is missing" );
    }

    for ( auto curIter : _mapOfConvertibleMaterialProperties ) {
        std::string nameOfProperty = curIter.second.getName();
        if ( curIter.second.hasValue() ) {
            nameOfProperty.resize( 16, ' ' );
            ( *char16Values )[position] = nameOfProperty.c_str();
            ( *doubleValues )[position] = curIter.second.getValue();
            ++position;
        }

        if ( curIter.second.isMandatory() && !curIter.second.hasValue() )
            throw std::runtime_error( "Mandatory material property " + nameOfProperty +
                                      " is missing" );
    }
    doubleValues->setSize( position );

    for ( auto curIter : _mapOfComplexMaterialProperties ) {
        std::string nameOfProperty = curIter.second.getName();
        if ( curIter.second.hasValue() ) {
            nameOfProperty.resize( 16, ' ' );
            ( *char16Values )[position] = nameOfProperty.c_str();
            ( *complexValues )[position] = curIter.second.getValue();
            ++position;
            ++pos3;
            if ( bOrdr )
                throw std::runtime_error( "ORDRE_PARAM with complex values not allowed" );
        }

        if ( curIter.second.isMandatory() && !curIter.second.hasValue() )
            throw std::runtime_error( "Mandatory material property " + nameOfProperty +
                                      " is missing" );
    }
    if ( pos3 != 0 )
        complexValues->setSize( position );
    else
        complexValues->setSize( pos3 );

    for ( auto curIter : _mapOfTableMaterialProperties ) {
        std::string nameOfProperty = curIter.second.getName();
        if ( curIter.second.hasValue() ) {
            nameOfProperty.resize( 16, ' ' );
            ( *char16Values )[position] = nameOfProperty.c_str();
            ( *char16Values )[position2] = curIter.second.getValue()->getName();
            ++position;
            ++position2;
        }

        if ( curIter.second.isMandatory() && !curIter.second.hasValue() )
            throw std::runtime_error( "Mandatory material property " + nameOfProperty +
                                      " is missing" );
    }

    for ( int i = 0; i < int(_vectKW.size()) / 2; ++i ) {
        auto nameOfProperty2 = _vectKW[2 * i];
        auto nameOfProperty = _vectKW[2 * i + 1];
        auto curIter1 = _mapOfFunctionMaterialProperties.find( nameOfProperty2 );
        if ( curIter1 == _mapOfFunctionMaterialProperties.end() )
            continue;

        auto curIter = *curIter1;
        if ( curIter.second.hasValue() ) {
            if ( bOrdr ) {
                const int curPos = 2 + 3 * nbOfMaterialProperties + i;
                ( *kOrd )[curPos] = position + 1;
                const int curPos2 = 2 + i;
                ( *kOrd )[curPos2] = mapTmp[nameOfProperty];
            }
            nameOfProperty.resize( 16, ' ' );
            ( *char16Values )[position] = nameOfProperty.c_str();
            ( *char16Values )[position2] = curIter.second.getValue()->getName();
            ++position;
            ++position2;
        }

        if ( curIter.second.isMandatory() && !curIter.second.hasValue() )
            throw std::runtime_error( "Mandatory material property " + nameOfProperty +
                                      " is missing" );
    }

    int position3 = 0;
    for ( auto curIter : _mapOfVectorRealMaterialProperties ) {
        std::string nameOfProperty = curIter.second.getName();
        if ( curIter.second.hasValue() ) {
            nameOfProperty.resize( 16, ' ' );
            ( *char16Values )[position] = nameOfProperty.c_str();
            ( *char16Values )[position2] = userReals[position3]->getName();

            auto values = curIter.second.getValue();
            userReals[position3]->allocate( values.size() );
            ( *userReals[position3] ) = values;
            ++position;
            ++position2;
            ++position3;
        }

        if ( curIter.second.isMandatory() && !curIter.second.hasValue() )
            throw std::runtime_error( "Mandatory material property " + nameOfProperty +
                                      " is missing" );
    }

    if ( _mapOfVectorFunctionMaterialProperties.size() > 1 )
        throw std::runtime_error( "Unconsistent size" );
    for ( auto curIter : _mapOfVectorFunctionMaterialProperties ) {
        std::string nameOfProperty = curIter.second.getName();
        if ( curIter.second.hasValue() ) {
            nameOfProperty.resize( 16, ' ' );
            ( *char16Values )[position] = nameOfProperty.c_str();
            ( *char16Values )[position2] = userFunctions[0]->getName();

            auto values = curIter.second.getValue();
            userFunctions[0]->allocate( values.size() );
            for ( int i = 0; i < int(values.size()); ++i )
                ( *userFunctions[0] )[i] = values[i]->getName();
            ++position;
            ++position2;
        }

        if ( curIter.second.isMandatory() && !curIter.second.hasValue() )
            throw std::runtime_error( "Mandatory material property " + nameOfProperty +
                                      " is missing" );
    }
    char16Values->setSize( position2 );

    return true;
};

bool GenericMaterialProperty::computeTractionFunction( FunctionPtr &doubleValues ) const
{
    return true;
};


int GenericMaterialProperty::getNumberOfPropertiesWithValue() const {
    int toReturn = 0;
    for ( auto curIter : _mapOfRealMaterialProperties )
        if ( curIter.second.hasValue() )
            ++toReturn;

    for ( auto curIter : _mapOfComplexMaterialProperties )
        if ( curIter.second.hasValue() )
            ++toReturn;

    for ( auto curIter : _mapOfStringMaterialProperties )
        if ( curIter.second.hasValue() )
            ++toReturn;

    for ( auto curIter : _mapOfTableMaterialProperties )
        if ( curIter.second.hasValue() )
            ++toReturn;

    for ( auto curIter : _mapOfFunctionMaterialProperties )
        if ( curIter.second.hasValue() )
            ++toReturn;

    for ( auto curIter : _mapOfVectorRealMaterialProperties )
        if ( curIter.second.hasValue() )
            ++toReturn;

    for ( auto curIter : _mapOfVectorFunctionMaterialProperties )
        if ( curIter.second.hasValue() )
            ++toReturn;

    for ( auto curIter : _mapOfConvertibleMaterialProperties )
        if ( curIter.second.hasValue() )
            ++toReturn;

    return toReturn;
};



bool GenericMaterialProperty::setValue( std::string nameOfProperty, ASTERDOUBLE value ) {
        // Recherche de la propriete materielle
        mapStrEMPDIterator curIter = _mapOfRealMaterialProperties.find( nameOfProperty );
        if ( curIter == _mapOfRealMaterialProperties.end() )
            return false;
        // Ajout de la valeur
        ( *curIter ).second.setValue( value );
        return true;
};

bool GenericMaterialProperty::setValue( std::string nameOfProperty, ASTERCOMPLEX value ) {
        // Recherche de la propriete materielle
        mapStrEMPCIterator curIter = _mapOfComplexMaterialProperties.find( nameOfProperty );
        if ( curIter == _mapOfComplexMaterialProperties.end() )
            return false;
        // Ajout de la valeur
        ( *curIter ).second.setValue( value );
        return true;
};


bool GenericMaterialProperty::setValue( std::string nameOfProperty, std::string value ) {
        // Recherche de la propriete materielle dans les Convertible
        const auto &curIter = _mapOfConvertibleMaterialProperties.find( nameOfProperty );
        if ( curIter != _mapOfConvertibleMaterialProperties.end() ) {
            ( *curIter ).second.setValue( value );
            return true;
        }

        // Recherche de la propriete materielle
        mapStrEMPSIterator curIter2 = _mapOfStringMaterialProperties.find( nameOfProperty );
        if ( curIter2 == _mapOfStringMaterialProperties.end() )
            return false;
        // Ajout de la valeur
        ( *curIter2 ).second.setValue( value );
        return true;
};


bool GenericMaterialProperty::setValue( std::string nameOfProperty, FunctionPtr value ) {
        // Recherche de la propriete materielle
        mapStrEMPFIterator curIter = _mapOfFunctionMaterialProperties.find( nameOfProperty );
        if ( curIter == _mapOfFunctionMaterialProperties.end() )
            return false;
        // Ajout de la valeur
        _mapOfFunctionMaterialProperties[nameOfProperty].setValue( value );
        return true;
};


bool GenericMaterialProperty::setValue( std::string nameOfProperty, TablePtr value ) {
        // Recherche de la propriete materielle
        mapStrEMPTIterator curIter = _mapOfTableMaterialProperties.find( nameOfProperty );
        if ( curIter == _mapOfTableMaterialProperties.end() )
            return false;
        // Ajout de la valeur
        ( *curIter ).second.setValue( value );
        return true;
};


bool GenericMaterialProperty::setValue( std::string nameOfProperty, Function2DPtr value ) {
        // Recherche de la propriete materielle
        mapStrEMPFIterator curIter = _mapOfFunctionMaterialProperties.find( nameOfProperty );
        if ( curIter == _mapOfFunctionMaterialProperties.end() )
            return false;
        // Ajout de la valeur
        _mapOfFunctionMaterialProperties[nameOfProperty].setValue( value );
        return true;
};

bool GenericMaterialProperty::setValue( std::string nameOfProperty, FormulaPtr value ) {
        // Recherche de la propriete materielle
        mapStrEMPFIterator curIter = _mapOfFunctionMaterialProperties.find( nameOfProperty );
        if ( curIter == _mapOfFunctionMaterialProperties.end() )
            return false;
        // Ajout de la valeur
        _mapOfFunctionMaterialProperties[nameOfProperty].setValue( value );
        return true;
};

bool GenericMaterialProperty::setValue( std::string nameOfProperty, VectorReal value ) {
        // Recherche de la propriete materielle
        auto curIter = _mapOfVectorRealMaterialProperties.find( nameOfProperty );
        if ( curIter == _mapOfVectorRealMaterialProperties.end() )
            return false;
        // Ajout de la valeur
        ( *curIter ).second._value = value;
        ( *curIter ).second._existsValue = true;
        return true;
};

bool GenericMaterialProperty::setValue( std::string nameOfProperty, VectorFunction value ) {
        // Recherche de la propriete materielle
        auto curIter = _mapOfVectorFunctionMaterialProperties.find( nameOfProperty );
        if ( curIter == _mapOfVectorFunctionMaterialProperties.end() )
            return false;
        // Ajout de la valeur
        ( *curIter ).second._value = value;
        ( *curIter ).second._existsValue = true;
        return true;
};
