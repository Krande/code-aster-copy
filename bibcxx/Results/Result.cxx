/**
 * @file Result.cxx
 * @brief Implementation de Result
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

#include "Results/Result.h"

#include "aster_fort_ds.h"
#include "aster_fort_jeveux.h"
#include "aster_fort_superv.h"
#include "aster_fort_utils.h"

#include "PythonBindings/LogicalUnitManager.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Tools.h"

std::pair< ASTERINTEGER, std::string > Result::_getNewFieldName( const std::string &name,
                                                                 const ASTERINTEGER &rank ) const {
    auto trim_name = trim( name );
    ASTERINTEGER retour;
    retour = 0;
    std::string null( " " );
    std::string returnName( 19, ' ' );
    CALLO_RSEXCH( null, getName(), trim_name, &rank, returnName, &retour );

    return std::make_pair( retour, returnName );
};

template < typename T >
void Result::_setFieldBase(
    const std::string &name, const ASTERINTEGER &rank, std::shared_ptr< T > field,
    std::map< std::string, std::map< ASTERINTEGER, std::shared_ptr< T > > > &dict ) {

    CALL_JEMARQ();

    if ( !field )
        raiseAsterError( "ValueError: field is empty" );

    auto trim_name = trim( name );
    auto rschex = _getNewFieldName( trim_name, rank );
    AS_ASSERT( rschex.first == 0 || rschex.first == 100 || rschex.first == 110 );

    if ( rschex.first == 110 ) {
        resize( std::max( (ASTERINTEGER)1, 2 * getNumberOfRanks() ) );
        rschex = _getNewFieldName( trim_name, rank );
    };

    CALLO_RSNOCH_FORWARD( getName(), name, &rank );
    std::string internalName( rschex.second.c_str(), 19 );

    if ( dict.count( trim_name ) == 0 ) {
        dict[trim_name] = std::map< ASTERINTEGER, std::shared_ptr< T > >();
    }

    // if field arlreday exist, destroy it befor to create new one
    if ( dict[trim_name].count( rank ) > 0 ) {
        dict[trim_name][rank] = nullptr;
    }

    auto result = std::make_shared< T >( internalName, *field );
    dict[trim_name][rank] = result;

    CALL_JEDEMA();
};

void Result::_checkMesh( const BaseMeshPtr mesh ) const {
    if ( !mesh )
        raiseAsterError( "ValueError: Mesh is empty" );

    if ( _mesh ) {
        if ( _mesh != mesh )
            raiseAsterError( "Incompatible meshes" );
    }
}

void Result::setMesh( const BaseMeshPtr &mesh ) {
    _checkMesh( mesh );
    _mesh = mesh;
};

void Result::setElementaryCharacteristics( const ElementaryCharacteristicsPtr &cara,
                                           ASTERINTEGER rank ) {
    if ( !cara )
        raiseAsterError( "ValueError: ElementaryCharacteristics is empty" );

    _mapElemCara[rank] = cara;
    std::string type( "CARAELEM" );
    CALLO_RSADPA_ZK8_WRAP( getName(), &rank, cara->getName(), type );
    setMesh( cara->getMesh() );
};

void Result::setListOfLoads( const ListOfLoadsPtr &load, ASTERINTEGER rank ) {
    if ( !load )
        raiseAsterError( "ValueError: Load is empty" );

    _mapLoads[rank] = load;
    std::string type( "EXCIT" );
    std::string cel( "E" );
    CALLO_RSADPA_ZK24_WRAP( &rank, getName(), load->getName(), type, cel );
};

void Result::setMaterialField( const MaterialFieldPtr &mater, ASTERINTEGER rank ) {
    if ( !mater )
        raiseAsterError( "ValueError: MaterialField is empty" );

    _mapMaterial[rank] = mater;
    std::string type( "CHAMPMAT" );
    CALLO_RSADPA_ZK8_WRAP( getName(), &rank, mater->getName(), type );
    setMesh( mater->getMesh() );
};

void Result::setModel( const ModelPtr &model, ASTERINTEGER rank ) {
    if ( !model )
        raiseAsterError( "ValueError: Model is empty" );

    _mapModel[rank] = model;
    std::string type( "MODELE" );
    CALLO_RSADPA_ZK8_WRAP( getName(), &rank, model->getName(), type );
    const auto fed = model->getFiniteElementDescriptor();
    _fieldBuidler.addFiniteElementDescriptor( fed );
    setMesh( model->getMesh() );
};

void Result::setParameterValue( std::string name, ASTERDOUBLE value, ASTERINTEGER rank ) {
    CALLO_RSADPA_ZR_WRAP( getName(), &rank, &value, name );
};

ASTERDOUBLE Result::getTimeValue( ASTERINTEGER rank ) {

    ASTERINTEGER nb_ranks = getNumberOfRanks();

    AS_ASSERT( rank <= nb_ranks );

    _rspr->updateValuePointer();

    for ( const auto &[i, item] : *_calculationParameter ) {
        item->updateValuePointer();
        auto typevar = trim( ( *item )[3].toString() );

        if ( typevar == "ACCES" ) {
            auto var_name = trim( _accessVariables->getStringFromIndex( i ) );
            if ( var_name == "INST" ) {
                auto nosuff = trim( ( *item )[0].toString() );
                auto ivar = std::stoi( trim( ( *item )[1].toString() ) );
                auto nmax = std::stoi( trim( ( *item )[2].toString() ) );

                AS_ASSERT( nosuff == ".RSPR" )

                auto index = nmax * rank + ivar - 1;

                return ( *_rspr )[index];
            }
        }
    }

    AS_ASSERT( false );
    return 0.0;
};

void Result::allocate( ASTERINTEGER nbRanks ) {

    std::string base( JeveuxMemoryTypesNames[Permanent] );
    ASTERINTEGER nbordr = nbRanks;
    CALLO_RSCRSD( base, getName(), getType(), &nbordr );

    AS_ASSERT( _calculationParameter->build( true ) );
    AS_ASSERT( _namesOfFields->build( true ) );
};

void Result::setElementaryCharacteristics( const ElementaryCharacteristicsPtr &cara ) {
    ASTERINTEGER nbRanks = getNumberOfRanks();
    if ( nbRanks > 0 ) {
        _serialNumber->updateValuePointer();
        for ( ASTERINTEGER rank = 0; rank < nbRanks; ++rank ) {
            const ASTERINTEGER iordr = ( *_serialNumber )[rank];
            if ( _mapElemCara.find( iordr ) == _mapElemCara.end() )
                setElementaryCharacteristics( cara, iordr );
        }
    }
};

void Result::setMaterialField( const MaterialFieldPtr &mater ) {
    ASTERINTEGER nbRanks = getNumberOfRanks();
    if ( nbRanks > 0 ) {
        _serialNumber->updateValuePointer();
        for ( ASTERINTEGER rank = 0; rank < nbRanks; ++rank ) {
            const ASTERINTEGER iordr = ( *_serialNumber )[rank];
            if ( _mapMaterial.find( iordr ) == _mapMaterial.end() )
                setMaterialField( mater, iordr );
        }
    }
};

void Result::setModel( const ModelPtr &model ) {
    ASTERINTEGER nbRanks = getNumberOfRanks();
    if ( nbRanks > 0 ) {
        _serialNumber->updateValuePointer();
        for ( ASTERINTEGER rank = 0; rank < nbRanks; ++rank ) {
            const ASTERINTEGER iordr = ( *_serialNumber )[rank];
            if ( _mapModel.find( iordr ) == _mapModel.end() )
                setModel( model, iordr );
        }
    }
};

std::vector< ElementaryCharacteristicsPtr > Result::getAllElementaryCharacteristics() const {
    return unique( _mapElemCara );
};

ElementaryCharacteristicsPtr Result::getElementaryCharacteristics() const {
    const auto cara = getAllElementaryCharacteristics();
    AS_ASSERT( cara.size() <= 1 );

    if ( cara.size() == 1 )
        return cara[0];

    return ElementaryCharacteristicsPtr( nullptr );
};

ElementaryCharacteristicsPtr Result::getElementaryCharacteristics( ASTERINTEGER rank ) const {
    return _mapElemCara.at( rank );
};

bool Result::hasElementaryCharacteristics( ASTERINTEGER rank ) const {
    return _mapElemCara.count( rank ) > 0;
};

bool Result::hasElementaryCharacteristics() const { return !_mapElemCara.empty(); };

bool Result::hasListOfLoads() const { return !_mapLoads.empty(); };

bool Result::hasListOfLoads( const ASTERINTEGER &rank ) const {
    return _mapLoads.count( rank ) > 0;
};

ListOfLoadsPtr Result::getListOfLoads( ASTERINTEGER rank ) const { return _mapLoads.at( rank ); };

std::vector< MaterialFieldPtr > Result::getMaterialFields() const {
    return unique( _mapMaterial );
};

MaterialFieldPtr Result::getMaterialField() const {
    const auto mate = getMaterialFields();
    AS_ASSERT( mate.size() <= 1 );

    if ( mate.size() == 1 )
        return mate[0];

    return MaterialFieldPtr( nullptr );
};

MaterialFieldPtr Result::getMaterialField( ASTERINTEGER rank ) const {
    return _mapMaterial.at( rank );
};

BaseMeshPtr Result::getMesh() const {
    if ( _mesh != nullptr )
        return _mesh;
    const auto model = getModel();
    if ( model != nullptr )
        return model->getMesh();
    return BaseMeshPtr( nullptr );
};

bool Result::hasMultipleModel() const {
    std::string name( "" );
    for ( const auto &curIter : _mapModel ) {
        if ( name == "" ) {
            name = curIter.second->getName();
        }
        if ( name != curIter.second->getName() )
            return true;
    }
    return false;
}

bool Result::hasModel( const ASTERINTEGER &rank ) const { return _mapModel.count( rank ) > 0; }

bool Result::hasMaterialField( const ASTERINTEGER &rank ) const {
    return _mapMaterial.count( rank ) > 0;
}

std::vector< ModelPtr > Result::getModels() const { return unique( _mapModel ); };

ModelPtr Result::getModel() const {
    if ( hasMultipleModel() ) {
        raiseAsterError( "Error: multiple models" );
    }

    const auto models = getModels();
    AS_ASSERT( models.size() <= 1 );

    if ( models.size() == 1 )
        return models[0];

    return ModelPtr( nullptr );
};

ModelPtr Result::getModel( ASTERINTEGER rank ) const { return _mapModel.at( rank ); };

ASTERINTEGER Result::getNumberOfRanks() const { return _serialNumber->size(); };

VectorLong Result::getRanks() const { return _serialNumber->toVector(); };

FieldOnCellsRealPtr Result::getFieldOnCellsReal( const std::string name,
                                                 const ASTERINTEGER rank ) const {
    return _dictOfMapOfFieldOnCellsReal.at( name ).at( rank );
};

FieldOnCellsComplexPtr Result::getFieldOnCellsComplex( const std::string name,
                                                       const ASTERINTEGER rank ) const {
    return _dictOfMapOfFieldOnCellsComplex.at( name ).at( rank );
};

FieldOnCellsLongPtr Result::getFieldOnCellsLong( const std::string name,
                                                 const ASTERINTEGER rank ) const {
    return _dictOfMapOfFieldOnCellsLong.at( name ).at( rank );
};

ConstantFieldOnCellsChar16Ptr
Result::getConstantFieldOnCellsChar16( const std::string name, const ASTERINTEGER rank ) const {
    return _dictOfMapOfConstantFieldOnCellsChar16.at( name ).at( rank );
};

ConstantFieldOnCellsRealPtr Result::getConstantFieldOnCellsReal( const std::string name,
                                                                 const ASTERINTEGER rank ) const {
    return _dictOfMapOfConstantFieldOnCellsReal.at( name ).at( rank );
};

py::dict Result::getAccessParameters() const {

    py::dict returnDict;
    std::string var_name, str_val, typevar, nosuff;
    ASTERINTEGER ivar, nmax, index;

    CALL_JEMARQ();
    _serialNumber->updateValuePointer();
    if ( _rspr->exists() )
        _rspr->updateValuePointer();
    if ( _rspi->exists() )
        _rspi->updateValuePointer();
    if ( _rsp8->exists() )
        _rsp8->updateValuePointer();
    if ( _rs16->exists() )
        _rs16->updateValuePointer();
    if ( _rs24->exists() )
        _rs24->updateValuePointer();

    AS_ASSERT( _calculationParameter->build() );

    ASTERINTEGER nb_ranks = getNumberOfRanks();

    var_name = "NUME_ORDRE";
    py::list listValues;
    for ( ASTERINTEGER j = 0; j < nb_ranks; ++j ) {
        listValues.append( ( *_serialNumber )[j] );
    }
    returnDict[var_name.c_str()] = listValues;

    auto items = _calculationParameter->getObjects();
    for ( auto& item : items ) {
        item->updateValuePointer();
        typevar = trim( (*item)[3].toString() );

        if ( typevar == "ACCES" ) {
            var_name = trim( _accessVariables->getStringFromIndex( item->getIndex() ) );
            nosuff = trim( ( *item )[0].toString() );
            ivar = std::stoi( trim( ( *item )[1].toString() ) );
            nmax = std::stoi( trim( ( *item )[2].toString() ) );

            py::list listV;

            if ( nosuff == ".RSPI" ) {
                for ( ASTERINTEGER j = 0; j < nb_ranks; ++j ) {
                    index = nmax * ( j ) + ivar - 1;
                    listV.append( ( *_rspi )[index] );
                }
            }

            else if ( nosuff == ".RSPR" ) {
                for ( ASTERINTEGER j = 0; j < nb_ranks; ++j ) {
                    index = nmax * ( j ) + ivar - 1;
                    listV.append( ( *_rspr )[index] );
                }
            }

            else {
                for ( ASTERINTEGER j = 0; j < nb_ranks; ++j ) {
                    index = nmax * ( j ) + ivar - 1;
                    if ( nosuff == ".RSP8" ) {
                        str_val = trim( ( ( *_rsp8 )[index] ).toString() );
                    } else if ( nosuff == ".RS16" ) {
                        str_val = trim( ( ( *_rs16 )[index] ).toString() );
                    } else if ( nosuff == ".RS24" ) {
                        str_val = trim( ( ( *_rs24 )[index] ).toString() );
                    } else {
                        AS_ASSERT( false );
                    }

                    if ( str_val.length() == 0 ) {
                        listV.append( py::none() );
                    } else {
                        listV.append( str_val );
                    }
                }
            }
            returnDict[var_name.c_str()] = listV;
        }
    }
    CALL_JEDEMA();
    return returnDict;
}

VectorString Result::getFieldsOnNodesRealNames() const {
    VectorString names;
    names.reserve( _dictOfMapOfFieldOnNodesReal.size() );

    for ( auto &it : _dictOfMapOfFieldOnNodesReal ) {
        std::string name = it.first;
        names.push_back( trim( name ) );
    }
    return names;
};

VectorString Result::getFieldsOnNodesComplexNames() const {
    VectorString names;
    names.reserve( _dictOfMapOfFieldOnNodesComplex.size() );

    for ( auto &it : _dictOfMapOfFieldOnNodesComplex ) {
        std::string name = it.first;
        names.push_back( trim( name ) );
    }
    return names;
};

VectorString Result::getFieldsOnCellsRealNames() const {
    VectorString names;
    names.reserve( _dictOfMapOfFieldOnCellsReal.size() );

    for ( auto &it : _dictOfMapOfFieldOnCellsReal ) {
        std::string name = it.first;
        names.push_back( trim( name ) );
    }
    return names;
};

VectorString Result::getFieldsOnCellsComplexNames() const {
    VectorString names;
    names.reserve( _dictOfMapOfFieldOnCellsComplex.size() );

    for ( auto &it : _dictOfMapOfFieldOnCellsComplex ) {
        std::string name = it.first;
        names.push_back( trim( name ) );
    }
    return names;
};

VectorString Result::getFieldsOnCellsLongNames() const {
    VectorString names;
    names.reserve( _dictOfMapOfFieldOnCellsLong.size() );

    for ( auto &it : _dictOfMapOfFieldOnCellsLong ) {
        std::string name = it.first;
        names.push_back( trim( name ) );
    }
    return names;
};

VectorString Result::getConstantFieldsOnCellsChar16Names() const {
    VectorString names;
    names.reserve( _dictOfMapOfConstantFieldOnCellsChar16.size() );

    for ( auto &it : _dictOfMapOfConstantFieldOnCellsChar16 ) {
        std::string name = it.first;
        names.push_back( trim( name ) );
    }
    return names;
};

VectorString Result::getConstantFieldsOnCellsRealNames() const {
    VectorString names;
    names.reserve( _dictOfMapOfConstantFieldOnCellsReal.size() );

    for ( auto &it : _dictOfMapOfConstantFieldOnCellsReal ) {
        std::string name = it.first;
        names.push_back( trim( name ) );
    }
    return names;
};

FieldOnNodesRealPtr Result::getFieldOnNodesReal( const std::string name,
                                                 const ASTERINTEGER rank ) const {

    return _dictOfMapOfFieldOnNodesReal.at( name ).at( rank );
};

FieldOnNodesComplexPtr Result::getFieldOnNodesComplex( const std::string name,
                                                       const ASTERINTEGER rank ) const {

    return _dictOfMapOfFieldOnNodesComplex.at( name ).at( rank );
};

void Result::setField( const FieldOnNodesRealPtr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    _setFieldBase( name, rank, field, _dictOfMapOfFieldOnNodesReal );
    _fieldBuidler.addFieldOnNodesDescription( field->getDescription() );
};

void Result::setField( const FieldOnNodesComplexPtr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    _setFieldBase( name, rank, field, _dictOfMapOfFieldOnNodesComplex );
    _fieldBuidler.addFieldOnNodesDescription( field->getDescription() );
};

void Result::setField( const FieldOnCellsRealPtr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    _setFieldBase( name, rank, field, _dictOfMapOfFieldOnCellsReal );
};

void Result::setField( const FieldOnCellsComplexPtr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    _setFieldBase( name, rank, field, _dictOfMapOfFieldOnCellsComplex );
};

void Result::setField( const FieldOnCellsLongPtr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    _setFieldBase( name, rank, field, _dictOfMapOfFieldOnCellsLong );
};

void Result::setField( const ConstantFieldOnCellsChar16Ptr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    _setFieldBase( name, rank, field, _dictOfMapOfConstantFieldOnCellsChar16 );
};

void Result::setField( const ConstantFieldOnCellsRealPtr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    _setFieldBase( name, rank, field, _dictOfMapOfConstantFieldOnCellsReal );
};

VectorString Result::getFieldsNames() const {
    VectorString vField;
    for ( auto curIter : _dictOfMapOfFieldOnNodesReal ) {
        vField.push_back( curIter.first );
    }
    for ( auto curIter : _dictOfMapOfFieldOnNodesComplex ) {
        vField.push_back( curIter.first );
    }
    for ( auto curIter : _dictOfMapOfFieldOnCellsReal ) {
        vField.push_back( curIter.first );
    }
    for ( auto curIter : _dictOfMapOfFieldOnCellsComplex ) {
        vField.push_back( curIter.first );
    }
    for ( auto curIter : _dictOfMapOfFieldOnCellsLong ) {
        vField.push_back( curIter.first );
    }
    for ( auto curIter : _dictOfMapOfConstantFieldOnCellsChar16 ) {
        vField.push_back( curIter.first );
    }
    for ( auto curIter : _dictOfMapOfConstantFieldOnCellsReal ) {
        vField.push_back( curIter.first );
    }
    return vField;
}

void Result::printListOfFields() const {
    auto vField = getFieldsNames();
    std::cout << "Content of DataStructure : ";
    for ( auto field : vField ) {
        std::cout << field << " - ";
    }
    std::cout << std::endl;
}

void Result::printInfo() const {
    ASTERINTEGER umess( 6 );
    CALLO_RSINFO( getName(), &umess );
}

void Result::printMedFile( const std::string fileName, std::string medName, bool local ) const {
    LogicalUnitFile a( fileName, Binary, New );
    ASTERINTEGER retour = a.getLogicalUnit();
    CommandSyntax cmdSt( "IMPR_RESU" );

    SyntaxMapContainer dict;
    dict.container["FORMAT"] = "MED";
    dict.container["UNITE"] = retour;

    if ( getMesh()->isParallel() ) {
        dict.container["PROC0"] = "NON";
        if ( !local )
            dict.container["FICHIER_UNIQUE"] = "OUI";
    } else
        dict.container["PROC0"] = "OUI";

    ListSyntaxMapContainer listeResu;
    SyntaxMapContainer dict2;
    dict2.container["RESULTAT"] = getName();
    dict2.container["TOUT_ORDRE"] = "OUI";
    if ( !medName.empty() )
        dict2.container["NOM_RESU_MED"] = medName.substr( 0, 8 );
    listeResu.push_back( dict2 );
    dict.container["RESU"] = listeResu;

    cmdSt.define( dict );

    ASTERINTEGER op = 39;
    CALL_EXECOP( &op );
};

bool Result::build( const std::vector< FiniteElementDescriptorPtr > feds,
                    const std::vector< FieldOnNodesDescriptionPtr > fnds ) {
    CALL_JEMARQ();
    _serialNumber->updateValuePointer();

    AS_ASSERT( _calculationParameter->build( true ) );
    AS_ASSERT( _namesOfFields->build( true ) );

    const auto nbRanks = getNumberOfRanks();

    for ( auto &fed : feds ) {
        _fieldBuidler.addFiniteElementDescriptor( fed );
    }

    for ( auto &fnd : fnds ) {
        _fieldBuidler.addFieldOnNodesDescription( fnd );
    }

    ASTERINTEGER cmpt = 1;
    for ( const auto &[key, obj] : _namesOfFields ) {
        obj->updateValuePointer();
        auto nomSymb = trim( _symbolicNamesOfFields->getStringFromIndex( cmpt ) );
        AS_ASSERT( nbRanks <= obj->size() );

        for ( ASTERINTEGER index = 0; index < nbRanks; ++index ) {
            std::string name( trim( ( *obj )[index].toString() ) );
            if ( name != "" ) {
                const ASTERINTEGER rank = ( *_serialNumber )[index];
                CALL_JEMARQ();
                std::string questi( "TYPE_CHAMP" );
                const std::string typeco( "CHAMP" );
                ASTERINTEGER repi = 0, ier = 0;
                JeveuxChar32 repk( " " );
                const std::string arret( "C" );

                CALLO_DISMOI( questi, name, typeco, &repi, repk, arret, &ier );
                const std::string resu( trim( repk.toString() ) );

                questi = "TYPE_SCA";
                repk = " ";
                CALLO_DISMOI( questi, name, typeco, &repi, repk, arret, &ier );
                const std::string scalaire( trim( repk.toString() ) );

                if ( resu == "NOEU" ) {
                    if ( scalaire == "R" ) {
                        if ( _dictOfMapOfFieldOnNodesReal.count( nomSymb ) == 0 ) {
                            _dictOfMapOfFieldOnNodesReal[nomSymb] = MapOfFieldOnNodesReal();
                        }

                        if ( _dictOfMapOfFieldOnNodesReal[nomSymb].count( rank ) == 0 ) {
                            AS_ASSERT( _mesh != nullptr );
                            FieldOnNodesRealPtr result =
                                _fieldBuidler.buildFieldOnNodes< ASTERDOUBLE >( name );
                            result->setMesh( _mesh );
                            _dictOfMapOfFieldOnNodesReal[nomSymb][rank] = result;
                        }
                    } else if ( scalaire == "C" ) {
                        if ( _dictOfMapOfFieldOnNodesComplex.count( nomSymb ) == 0 ) {
                            _dictOfMapOfFieldOnNodesComplex[nomSymb] = MapOfFieldOnNodesComplex();
                        }

                        if ( _dictOfMapOfFieldOnNodesComplex[nomSymb].count( rank ) == 0 ) {
                            AS_ASSERT( _mesh != nullptr );
                            FieldOnNodesComplexPtr result =
                                _fieldBuidler.buildFieldOnNodes< ASTERCOMPLEX >( name );
                            result->setMesh( _mesh );
                            _dictOfMapOfFieldOnNodesComplex[nomSymb][rank] = result;
                        }
                    } else {
                        AS_ABORT( "Type not supported: " + scalaire );
                    }

                } else if ( resu == "ELEM" || resu == "ELNO" || resu == "ELGA" ) {
                    if ( scalaire == "R" ) {
                        if ( _dictOfMapOfFieldOnCellsReal.count( nomSymb ) == 0 ) {
                            _dictOfMapOfFieldOnCellsReal[nomSymb] = MapOfFieldOnCellsReal();
                        }

                        if ( _dictOfMapOfFieldOnCellsReal[nomSymb].count( rank ) == 0 ) {
                            AS_ASSERT( _mesh != nullptr );
                            auto result =
                                _fieldBuidler.buildFieldOnCells< ASTERDOUBLE >( name, _mesh );
                            _dictOfMapOfFieldOnCellsReal[nomSymb][rank] = result;
                        }
                    } else if ( scalaire == "C" ) {
                        if ( _dictOfMapOfFieldOnCellsComplex.count( nomSymb ) == 0 ) {
                            _dictOfMapOfFieldOnCellsComplex[nomSymb] = MapOfFieldOnCellsComplex();
                        }

                        if ( _dictOfMapOfFieldOnCellsComplex[nomSymb].count( rank ) == 0 ) {
                            AS_ASSERT( _mesh != nullptr );
                            auto result =
                                _fieldBuidler.buildFieldOnCells< ASTERCOMPLEX >( name, _mesh );
                            _dictOfMapOfFieldOnCellsComplex[nomSymb][rank] = result;
                        }
                    } else if ( scalaire == "I" ) {
                        if ( _dictOfMapOfFieldOnCellsLong.count( nomSymb ) == 0 ) {
                            _dictOfMapOfFieldOnCellsLong[nomSymb] = MapOfFieldOnCellsLong();
                        }

                        if ( _dictOfMapOfFieldOnCellsLong[nomSymb].count( rank ) == 0 ) {
                            AS_ASSERT( _mesh != nullptr );
                            auto result =
                                _fieldBuidler.buildFieldOnCells< ASTERINTEGER >( name, _mesh );
                            _dictOfMapOfFieldOnCellsLong[nomSymb][rank] = result;
                        }
                    } else {
                        AS_ABORT( "Type not supported: " + scalaire );
                    }
                } else if ( resu == "CART" ) {
                    if ( scalaire == "K16" ) {
                        if ( _dictOfMapOfConstantFieldOnCellsChar16.count( nomSymb ) == 0 ) {
                            _dictOfMapOfConstantFieldOnCellsChar16[nomSymb] =
                                MapOfConstantFieldOnCellsChar16();
                        }

                        if ( _dictOfMapOfConstantFieldOnCellsChar16[nomSymb].count( rank ) == 0 ) {
                            AS_ASSERT( _mesh != nullptr );
                            auto result = _fieldBuidler.buildConstantFieldOnCells< JeveuxChar16 >(
                                name, _mesh );
                            _dictOfMapOfConstantFieldOnCellsChar16[nomSymb][rank] = result;
                        }
                    } else if ( scalaire == "R" ) {
                        if ( _dictOfMapOfConstantFieldOnCellsReal.count( nomSymb ) == 0 ) {
                            _dictOfMapOfConstantFieldOnCellsReal[nomSymb] =
                                MapOfConstantFieldOnCellsReal();
                        }

                        if ( _dictOfMapOfConstantFieldOnCellsReal[nomSymb].count( rank ) == 0 ) {
                            AS_ASSERT( _mesh != nullptr );
                            auto result = _fieldBuidler.buildConstantFieldOnCells< ASTERDOUBLE >(
                                name, _mesh );
                            _dictOfMapOfConstantFieldOnCellsReal[nomSymb][rank] = result;
                        }
                    } else {
                        AS_ABORT( "Type not supported: " + scalaire );
                    }
                } else {
                    std::cout << "Field not build : " << name << " (" << resu << ")" << std::endl;
                }
                CALL_JEDEMA();
            }
        }
        ++cmpt;
    }

    CALL_JEDEMA();
    return update_tables();
};

void Result::resize( ASTERINTEGER nbRanks ) {
    if ( getNumberOfRanks() == 0 ) {
        allocate( nbRanks );
    } else {
        CALLO_RSAGSD( getName(), &nbRanks );
    }
}

std::vector< FiniteElementDescriptorPtr > Result::getFiniteElementDescriptors() const {
    return _fieldBuidler.getFiniteElementDescriptors();
};

std::vector< FieldOnNodesDescriptionPtr > Result::getFieldOnNodesDescriptions() const {
    return _fieldBuidler.getFieldOnNodesDescriptions();
};
