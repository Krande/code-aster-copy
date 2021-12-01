/**
 * @file TableContainer.cxx
 * @brief Implementation de TableContainer
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

#include "DataFields/TableContainer.h"
#include <iostream>

/* person_in_charge: nicolas.sellenet at edf.fr */

void TableContainer::addObject( const std::string &a, ElementaryMatrixDisplacementRealPtr b ) {
    _mapEMDD[a] = b;
};

void TableContainer::addObject( const std::string &a, ElementaryMatrixTemperatureRealPtr b ) {
    _mapEMTD[a] = b;
};

void TableContainer::addObject( const std::string &a, ElementaryVectorDisplacementRealPtr b ) {
    _mapEVDD[a] = b;
};

void TableContainer::addObject( const std::string &a, ElementaryVectorTemperatureRealPtr b ) {
    _mapEVTD[a] = b;
};

void TableContainer::addObject( const std::string &a, FieldOnCellsRealPtr b ) {
    _mapFOED[a] = b;
};

void TableContainer::addObject( const std::string &a, FieldOnNodesRealPtr b ) {
    _mapFOND[a] = b;
};

void TableContainer::addObject( const std::string &a, FunctionPtr b ) { _mapF[a] = b; };

void TableContainer::addObject( const std::string &a, FunctionComplexPtr b ) {
    _mapFC[a] = b;
};

void TableContainer::addObject( const std::string &a, GeneralizedAssemblyMatrixRealPtr b ) {
    _mapGAMD[a] = b;
};

void TableContainer::addObject( const std::string &a, DataFieldPtr b ) { _mapGDF[a] = b; };

void TableContainer::addObject( const std::string &a, ModeResultPtr b ) { _mapMMC[a] = b; };

void TableContainer::addObject( const std::string &a, ConstantFieldOnCellsRealPtr b ) {
    _mapPCFOMD[a] = b;
};

void TableContainer::addObject( const std::string &a, Function2DPtr b ) { _mapS[a] = b; };

void TableContainer::addObject( const std::string &a, TablePtr b ) { _mapT[a] = b; };

ElementaryMatrixDisplacementRealPtr
TableContainer::getElementaryMatrixDisplacementReal( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapEMDD.find( aa );
    if ( curIter == _mapEMDD.end() )
        return ElementaryMatrixDisplacementRealPtr( nullptr );
    return curIter->second;
};

ElementaryMatrixTemperatureRealPtr
TableContainer::getElementaryMatrixTemperatureReal( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapEMTD.find( aa );
    if ( curIter == _mapEMTD.end() )
        return ElementaryMatrixTemperatureRealPtr( nullptr );
    return curIter->second;
};

ElementaryVectorDisplacementRealPtr
TableContainer::getElementaryVectorDisplacementReal( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapEVDD.find( aa );
    if ( curIter == _mapEVDD.end() )
        return ElementaryVectorDisplacementRealPtr( nullptr );
    return curIter->second;
};

ElementaryVectorTemperatureRealPtr
TableContainer::getElementaryVectorTemperatureReal( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapEVTD.find( aa );
    if ( curIter == _mapEVTD.end() )
        return ElementaryVectorTemperatureRealPtr( nullptr );
    return curIter->second;
};

FieldOnCellsRealPtr TableContainer::getFieldOnCellsReal( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapFOED.find( aa );
    if ( curIter == _mapFOED.end() )
        return FieldOnCellsRealPtr( nullptr );
    return curIter->second;
};

FieldOnNodesRealPtr TableContainer::getFieldOnNodesReal( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapFOND.find( aa );
    if ( curIter == _mapFOND.end() )
        return FieldOnNodesRealPtr( nullptr );
    return curIter->second;
};

FunctionPtr TableContainer::getFunction( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapF.find( aa );
    if ( curIter == _mapF.end() )
        return FunctionPtr( nullptr );
    return curIter->second;
};

FunctionComplexPtr TableContainer::getFunctionComplex( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapFC.find( aa );
    if ( curIter == _mapFC.end() )
        return FunctionComplexPtr( nullptr );
    return curIter->second;
};

GeneralizedAssemblyMatrixRealPtr
TableContainer::getGeneralizedAssemblyMatrix( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapGAMD.find( aa );
    if ( curIter == _mapGAMD.end() )
        return GeneralizedAssemblyMatrixRealPtr( nullptr );
    return curIter->second;
};

DataFieldPtr TableContainer::getDataField( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapGDF.find( aa );
    if ( curIter == _mapGDF.end() )
        return DataFieldPtr( nullptr );
    return curIter->second;
};

ModeResultPtr TableContainer::getModeResult( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapMMC.find( aa );
    if ( curIter == _mapMMC.end() )
        return ModeResultPtr( nullptr );
    return curIter->second;
};

ConstantFieldOnCellsRealPtr
TableContainer::getConstantFieldOnCellsReal( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapPCFOMD.find( aa );
    if ( curIter == _mapPCFOMD.end() )
        return ConstantFieldOnCellsRealPtr( nullptr );
    return curIter->second;
};

Function2DPtr TableContainer::getFunction2D( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapS.find( aa );
    if ( curIter == _mapS.end() )
        return Function2DPtr( nullptr );
    return curIter->second;
};

TablePtr TableContainer::getTable( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapT.find( aa );
    if ( curIter == _mapT.end() )
        return TablePtr( nullptr );
    return curIter->second;
};

bool TableContainer::build() {

    _parameterDescription->updateValuePointer();
    const int size = _parameterDescription->size() / 4;
    for ( int i = 0; i < size; ++i ) {
        const auto test = trim( ( *_parameterDescription )[i * 4].toString() );
        const auto name = trim( ( *_parameterDescription )[i * 4 + 2].toString() );
        const auto name2 = trim( ( *_parameterDescription )[i * 4 + 3].toString() );
        const auto type = trim( ( *_parameterDescription )[i * 4 + 1].toString() );
        if ( test == "NOM_OBJET" ) {
            if ( _objectName.isEmpty() )
                _objectName = JeveuxVectorChar16( name );
        } else if ( test == "TYPE_OBJET" ) {
            if ( _objectType.isEmpty() )
                _objectType = JeveuxVectorChar16( name );
        } else if ( test == "NOM_SD" ) {
            if ( type == "K8" ) {
                if ( _dsName1.isEmpty() )
                    _dsName1 = JeveuxVectorChar8( name );
            } else {
                if ( _dsName2.isEmpty() )
                    _dsName2 = JeveuxVectorChar24( name );
            }
        } else
            _others.push_back( JeveuxVectorLong( name ) );
        _vecOfSizes.push_back( JeveuxVectorLong( name2 ) );
    }
    const bool test1 = _dsName1.isEmpty();
    const bool test2 = _dsName2.isEmpty();
    const bool test3 = test1 && test2;
    if ( !test3 && !_objectType.isEmpty() && !_objectName.isEmpty() ) {
        int usedSize = 0;
        if ( !test1 ) {
            _dsName1->updateValuePointer();
            usedSize = _dsName1->size();
        } else {
            _dsName2->updateValuePointer();
            usedSize = _dsName2->size();
        }
        _objectType->updateValuePointer();
        _objectName->updateValuePointer();

        if ( usedSize != _objectType->size() )
            throw std::runtime_error( "Unconsistent sizes" );
        if ( usedSize != _objectName->size() )
            throw std::runtime_error( "Unconsistent sizes" );
        for ( int i = 0; i < usedSize; ++i ) {
            const auto type = trim( ( *_objectType )[i].toString() );
            std::string sdName( "" );
            if ( !test1 )
                sdName = trim( ( *_dsName1 )[i].toString() );
            else
                sdName = trim( ( *_dsName2 )[i].toString() );
            const auto name = trim( ( *_objectName )[i].toString() );
            if ( type.empty() && sdName.empty() ) {
                // pass
            } else if ( type == "MATR_ASSE_GENE_R" ) {
                if ( _mapGAMD[name] == nullptr )
                    _mapGAMD[name] = GeneralizedAssemblyMatrixRealPtr(
                        new GeneralizedAssemblyMatrixReal( sdName ) );
            } else if ( type == "MATR_ELEM_DEPL_R" ) {
                if ( _mapEMDD[name] == nullptr )
                    _mapEMDD[name] = ElementaryMatrixDisplacementRealPtr(
                        new ElementaryMatrixDisplacementReal( sdName ) );
            } else if ( type == "MATR_ELEM_TEMP_R" ) {
                if ( _mapEMTD[name] == nullptr )
                    _mapEMTD[name] = ElementaryMatrixTemperatureRealPtr(
                        new ElementaryMatrixTemperatureReal( sdName ) );
            } else if ( type == "VECT_ELEM_DEPL_R" ) {
                if ( _mapEVDD[name] == nullptr )
                    _mapEVDD[name] = ElementaryVectorDisplacementRealPtr(
                        new ElementaryVectorDisplacementReal( sdName ) );
            } else if ( type == "VECT_ELEM_TEMP_R" ) {
                if ( _mapEVTD[name] == nullptr )
                    _mapEVTD[name] = ElementaryVectorTemperatureRealPtr(
                        new ElementaryVectorTemperatureReal( sdName ) );
            } else if ( type == "CHAM_GD_SDASTER" ) {
                if ( _mapGDF[name] == nullptr )
                    _mapGDF[name] = DataFieldPtr( new DataField( sdName, "CHAM_GD" ) );
            } else if ( type == "CHAM_NO_SDASTER" ) {
                if ( _mapFOND[name] == nullptr )
                    _mapFOND[name] = FieldOnNodesRealPtr( new FieldOnNodesReal( sdName ) );
            }
            //             else if( type == "CARTE_SDASTER" )
            //                 _mapPCFOMD[name] = ConstantFieldOnCellsRealPtr
            //                                     ( new ConstantFieldOnCellsReal( sdName ) );
            else if ( type == "CHAM_ELEM" ) {
                if ( _mapFOED[name] == nullptr )
                    _mapFOED[name] = FieldOnCellsRealPtr( new FieldOnCellsReal( sdName ) );
            } else if ( type == "MODE_MECA" ) {
                if ( _mapMMC[name] == nullptr )
                    _mapMMC[name] = ModeResultPtr( new ModeResult( sdName ) );
            } else if ( type == "TABLE_SDASTER" ) {
                if ( _mapT[name] == nullptr )
                    _mapT[name] = TablePtr( new Table( sdName ) );
            } else if ( type == "FONCTION_SDASTER" ) {
                if ( _mapF[name] == nullptr )
                    _mapF[name] = FunctionPtr( new Function( sdName ) );
            } else if ( type == "FONCTION_C" ) {
                if ( _mapFC[name] == nullptr )
                    _mapFC[name] = FunctionComplexPtr( new FunctionComplex( sdName ) );
            } else if ( type == "NAPPE_SDASTER" ) {
                if ( _mapS[name] == nullptr )
                    _mapS[name] = Function2DPtr( new Function2D( sdName ) );
            } else
                throw std::runtime_error( "Type not implemented '" + type + "' for '" + name +
                                          "'" );
        }
    } else
        return false;
    return true;
};
