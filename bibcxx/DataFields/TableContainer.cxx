/**
 * @file TableContainer.cxx
 * @brief Implementation de TableContainer
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

void TableContainer::addObject( const std::string &a, FieldOnCellsRealPtr b ) { _mapFOED[a] = b; };

void TableContainer::addObject( const std::string &a, FieldOnNodesRealPtr b ) { _mapFOND[a] = b; };

void TableContainer::addObject( const std::string &a, MeshPtr b ) { _mapMesh[a] = b; };

#ifdef ASTER_HAVE_MPI
void TableContainer::addObject( const std::string &a, ParallelMeshPtr b ) { _mapPMesh[a] = b; };
#endif

void TableContainer::addObject( const std::string &a, FunctionPtr b ) { _mapF[a] = b; };

void TableContainer::addObject( const std::string &a, FunctionComplexPtr b ) { _mapFC[a] = b; };

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

MeshPtr TableContainer::getMesh( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapMesh.find( aa );
    if ( curIter == _mapMesh.end() )
        return MeshPtr( nullptr );
    return curIter->second;
};

#ifdef ASTER_HAVE_MPI
ParallelMeshPtr TableContainer::getParallelMesh( const std::string &a ) const {
    const auto aa = trim( a );
    const auto curIter = _mapPMesh.find( aa );
    if ( curIter == _mapPMesh.end() )
        return ParallelMeshPtr( nullptr );
    return curIter->second;
};
#endif

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
        const auto para = trim( ( *_parameterDescription )[i * 4].toString() );
        const auto type = trim( ( *_parameterDescription )[i * 4 + 1].toString() );
        const auto objev = trim( ( *_parameterDescription )[i * 4 + 2].toString() );

#ifdef ASTER_DEBUG_CXX
        std::cout << "DEBUG descr: para: " << para << " objev: " << objev << " type: " << type
                  << std::endl;
#endif
        if ( para == "NOM_OBJET" ) {
            if ( _objectName.isEmpty() )
                _objectName = JeveuxVectorChar16( objev );
        } else if ( para == "TYPE_OBJET" ) {
            if ( _objectType.isEmpty() )
                _objectType = JeveuxVectorChar16( objev );
        } else if ( para == "NOM_SD" ) {
            if ( type == "K8" ) {
                if ( _dsName8.isEmpty() )
                    _dsName8 = JeveuxVectorChar8( objev );
            } else {
                if ( _dsName24.isEmpty() )
                    _dsName24 = JeveuxVectorChar24( objev );
            }
        }
    }
    const bool usek8 = !_dsName8.isEmpty();
    const bool usek24 = !_dsName24.isEmpty();
    const bool hasDS = usek8 || usek24;
    AS_ASSERT( !( usek8 && usek24 ) );
    if ( hasDS && !_objectType.isEmpty() && !_objectName.isEmpty() ) {
        int usedSize = 0;
        if ( usek8 ) {
            _dsName8->updateValuePointer();
            usedSize = _dsName8->size();
        } else {
            _dsName24->updateValuePointer();
            usedSize = _dsName24->size();
        }
        _objectType->updateValuePointer();
        _objectName->updateValuePointer();

        if ( usedSize != _objectType->size() )
            throw std::runtime_error( "Unconsistent size for types" );
        if ( usedSize != _objectName->size() )
            throw std::runtime_error( "Unconsistent size for names" );
        for ( int i = 0; i < usedSize; ++i ) {
            std::string type = trim( ( *_objectType )[i].toString() );
            std::string dsName( "" );
            if ( usek8 )
                dsName = trim( ( *_dsName8 )[i].toString() );
            else
                dsName = trim( ( *_dsName24 )[i].toString() );
            const auto name = trim( ( *_objectName )[i].toString() );

#ifdef ASTER_DEBUG_CXX
            std::cout << "DEBUG: TableContainer index: " << i << " dsName: " << dsName
                      << " objName: " << name << " objType:" << type << std::endl;
#endif
            auto pos = type.find( "_SDASTER" );
            if ( pos ) {
                type = type.substr( 0, pos );
            }
            if ( type.empty() || dsName.empty() ) {
                // pass
            } else if ( type == "MATR_ASSE_GENE_R" ) {
                if ( _mapGAMD[name] == nullptr ) {
                    _mapGAMD[name] = std::make_shared< GeneralizedAssemblyMatrixReal >( dsName );
                }
            } else if ( type == "MATR_ELEM_DEPL_R" ) {
                if ( _mapEMDD[name] == nullptr ) {
                    _mapEMDD[name] = std::make_shared< ElementaryMatrixDisplacementReal >( dsName );
                }
            } else if ( type == "MATR_ELEM_TEMP_R" ) {
                if ( _mapEMTD[name] == nullptr ) {
                    _mapEMTD[name] = std::make_shared< ElementaryMatrixTemperatureReal >( dsName );
                }
            } else if ( type == "VECT_ELEM_DEPL_R" ) {
                if ( _mapEVDD[name] == nullptr ) {
                    _mapEVDD[name] = std::make_shared< ElementaryVectorDisplacementReal >( dsName );
                }
            } else if ( type == "VECT_ELEM_TEMP_R" ) {
                if ( _mapEVTD[name] == nullptr ) {
                    _mapEVTD[name] = std::make_shared< ElementaryVectorTemperatureReal >( dsName );
                }
            } else if ( type == "CHAM_GD" ) {
                if ( _mapGDF[name] == nullptr ) {
                    _mapGDF[name] = std::make_shared< DataField >( dsName, "CHAM_GD" );
                }
            } else if ( type == "CHAM_NO" ) {
                if ( _mapFOND[name] == nullptr ) {
                    _mapFOND[name] = std::make_shared< FieldOnNodesReal >( dsName );
                    _mapFOND[name]->build();
                }
                // } else if ( type == "CARTE" ) {
                //     _mapPCFOMD[name] = std::make_shared< ConstantFieldOnCellsReal >( dsName );
            } else if ( type == "CHAM_ELEM" ) {
                if ( _mapFOED[name] == nullptr ) {
                    _mapFOED[name] = std::make_shared< FieldOnCellsReal >( dsName );
                }
            } else if ( type == "MODE_MECA" ) {
                if ( _mapMMC[name] == nullptr ) {
                    _mapMMC[name] = std::make_shared< ModeResult >( dsName );
                }
            } else if ( type == "TABLE" ) {
                if ( _mapT[name] == nullptr ) {
                    _mapT[name] = std::make_shared< Table >( dsName );
                }
            } else if ( type == "MAILLAGE" ) {
                if ( _mapMesh[name] == nullptr ) {
                    _mapMesh[name] = std::make_shared< Mesh >( dsName );
                }
#ifdef ASTER_HAVE_MPI
            } else if ( type == "MAILLAGE_P" ) {
                if ( _mapPMesh[name] == nullptr ) {
                    _mapPMesh[name] = std::make_shared< ParallelMesh >( dsName );
                }
#endif
            } else if ( type == "FONCTION" ) {
                if ( _mapF[name] == nullptr ) {
                    _mapF[name] = std::make_shared< Function >( dsName );
                }
            } else if ( type == "FONCTION_C" ) {
                if ( _mapFC[name] == nullptr ) {
                    _mapFC[name] = std::make_shared< FunctionComplex >( dsName );
                }
            } else if ( type == "NAPPE" ) {
                if ( _mapS[name] == nullptr ) {
                    _mapS[name] = std::make_shared< Function2D >( dsName );
                }
            } else {
                throw std::runtime_error( "Unsupported type '" + type + "' for '" + name + "'" );
            }
        }
    } else
        return false;
    return true;
};
