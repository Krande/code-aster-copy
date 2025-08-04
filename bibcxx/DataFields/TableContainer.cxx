/**
 * @file TableContainer.cxx
 * @brief Implementation de TableContainer
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

#include "DataFields/TableContainer.h"

#include <iostream>

/* person_in_charge: nicolas.sellenet at edf.fr */

void TableContainer::addObject( const std::string &a, ElementaryMatrixDisplacementRealPtr b ) {
    _mapEMDD[a] = b;
};

void TableContainer::addObject( const std::string &a, ElementaryVectorDisplacementRealPtr b ) {
    _mapEVDD[a] = b;
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
    const auto aa = strip( a );
    const auto curIter = _mapEMDD.find( aa );
    if ( curIter == _mapEMDD.end() )
        return ElementaryMatrixDisplacementRealPtr( nullptr );
    return curIter->second;
};

ElementaryVectorDisplacementRealPtr
TableContainer::getElementaryVectorDisplacementReal( const std::string &a ) const {
    const auto aa = strip( a );
    const auto iter = _mapEVDD.find( aa );
    if ( iter != _mapEVDD.end() ) {
        return _mapEVDD.at( aa );
    }

    return ElementaryVectorDisplacementRealPtr( nullptr );
};

FieldOnCellsRealPtr TableContainer::getFieldOnCellsReal( const std::string &a ) const {
    const auto aa = strip( a );
    const auto curIter = _mapFOED.find( aa );
    if ( curIter == _mapFOED.end() )
        return FieldOnCellsRealPtr( nullptr );
    return curIter->second;
};

FieldOnNodesRealPtr TableContainer::getFieldOnNodesReal( const std::string &a ) const {
    const auto aa = strip( a );
    const auto curIter = _mapFOND.find( aa );
    if ( curIter == _mapFOND.end() )
        return FieldOnNodesRealPtr( nullptr );
    return curIter->second;
};

MeshPtr TableContainer::getMesh( const std::string &a ) const {
    const auto aa = strip( a );
    const auto curIter = _mapMesh.find( aa );
    if ( curIter == _mapMesh.end() )
        return MeshPtr( nullptr );
    return curIter->second;
};

#ifdef ASTER_HAVE_MPI
ParallelMeshPtr TableContainer::getParallelMesh( const std::string &a ) const {
    const auto aa = strip( a );
    const auto curIter = _mapPMesh.find( aa );
    if ( curIter == _mapPMesh.end() )
        return ParallelMeshPtr( nullptr );
    return curIter->second;
};
#endif

FunctionPtr TableContainer::getFunction( const std::string &a ) const {
    const auto aa = strip( a );
    const auto curIter = _mapF.find( aa );
    if ( curIter == _mapF.end() )
        return FunctionPtr( nullptr );
    return curIter->second;
};

FunctionComplexPtr TableContainer::getFunctionComplex( const std::string &a ) const {
    const auto aa = strip( a );
    const auto curIter = _mapFC.find( aa );
    if ( curIter == _mapFC.end() )
        return FunctionComplexPtr( nullptr );
    return curIter->second;
};

GeneralizedAssemblyMatrixRealPtr
TableContainer::getGeneralizedAssemblyMatrix( const std::string &a ) const {
    const auto aa = strip( a );
    const auto curIter = _mapGAMD.find( aa );
    if ( curIter == _mapGAMD.end() )
        return GeneralizedAssemblyMatrixRealPtr( nullptr );
    return curIter->second;
};

DataFieldPtr TableContainer::getDataField( const std::string &a ) const {
    const auto aa = strip( a );
    const auto curIter = _mapGDF.find( aa );
    if ( curIter == _mapGDF.end() )
        return DataFieldPtr( nullptr );
    return curIter->second;
};

ModeResultPtr TableContainer::getModeResult( const std::string &a ) const {
    const auto aa = strip( a );
    const auto curIter = _mapMMC.find( aa );
    if ( curIter == _mapMMC.end() )
        return ModeResultPtr( nullptr );
    return curIter->second;
};

ConstantFieldOnCellsRealPtr
TableContainer::getConstantFieldOnCellsReal( const std::string &a ) const {
    const auto aa = strip( a );
    const auto curIter = _mapPCFOMD.find( aa );
    if ( curIter == _mapPCFOMD.end() )
        return ConstantFieldOnCellsRealPtr( nullptr );
    return curIter->second;
};

Function2DPtr TableContainer::getFunction2D( const std::string &a ) const {
    const auto aa = strip( a );
    const auto curIter = _mapS.find( aa );
    if ( curIter == _mapS.end() )
        return Function2DPtr( nullptr );
    return curIter->second;
};

TablePtr TableContainer::getTable( const std::string &a ) const {
    const auto aa = strip( a );
    const auto curIter = _mapT.find( aa );
    if ( curIter == _mapT.end() )
        return TablePtr( nullptr );
    return curIter->second;
};

bool TableContainer::build() {
    Table::build();

    VectorString parameters = getParameters();
    for ( std::string parameter : { "NOM_SD", "NOM_OBJET", "TYPE_OBJET" } ) {
        if ( std::find( parameters.begin(), parameters.end(), parameter ) == parameters.end() )
            throw std::runtime_error( "missing parameter " + parameter + " in TableContainer" );
    }

    bool is_dsname_K24;
    int usedSize;
    std::string type = getColumnType( "NOM_SD" );
    if ( type == "K8" ) {
        _dsName = _columnChar8.at( "NOM_SD" );
        _dsName->updateValuePointer();
        usedSize = _dsName->size();
        is_dsname_K24 = false;
    } else if ( type == "K24" ) {
        _dsName24 = _columnChar24.at( "NOM_SD" );
        _dsName24->updateValuePointer();
        usedSize = _dsName24->size();
        is_dsname_K24 = true;
    } else
        AS_ASSERT( false )

    type = getColumnType( "NOM_OBJET" );
    AS_ASSERT( type == "K16" );
    _objectName = _columnChar16.at( "NOM_OBJET" );
    _objectName->updateValuePointer();
    if ( usedSize != _objectName->size() )
        throw std::runtime_error( "Unconsistent size for names" );

    bool is_typeobject_K24;
    type = getColumnType( "TYPE_OBJET" );
    if ( type == "K16" ) {
        _objectType = _columnChar16.at( "TYPE_OBJET" );
        _objectType->updateValuePointer();
        is_typeobject_K24 = false;
        if ( usedSize != _objectType->size() )
            throw std::runtime_error( "Unconsistent size for types" );
    } else if ( type == "K24" ) {
        _objectType24 = _columnChar24.at( "TYPE_OBJET" );
        _objectType24->updateValuePointer();
        is_typeobject_K24 = true;
        if ( usedSize != _objectType24->size() )
            throw std::runtime_error( "Unconsistent size for types" );
    } else
        AS_ASSERT( false );

    for ( int i = 0; i < usedSize; ++i ) {
        std::string type;
        if ( is_typeobject_K24 )
            type = strip( ( *_objectType24 )[i].toString() );
        else
            type = strip( ( *_objectType )[i].toString() );
        std::string dsName;
        if ( is_dsname_K24 )
            dsName = strip( ( *_dsName24 )[i].toString() );
        else
            dsName = strip( ( *_dsName )[i].toString() );
        std::string name = strip( ( *_objectName )[i].toString() );

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
        } else if ( type == "VECT_ELEM_DEPL_R" ) {
            const auto iter = _mapEVDD.find( name );
            if ( ( iter == _mapEVDD.end() ) || _mapEVDD.at( name ) == nullptr ) {
                _mapEVDD[name] =
                    std::make_shared< ElementaryVectorDisplacementReal >( dsName, nullptr );
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
                _mapT[name]->build();
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

    return true;
};
