/**
 * @file Result.cxx
 * @brief Implementation de Result
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

#include "aster_fort_ds.h"
#include "aster_fort_jeveux.h"
#include "aster_fort_superv.h"
#include "aster_fort_utils.h"

#include "PythonBindings/LogicalUnitManager.h"
#include "Results/Result.h"
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

void Result::_checkMesh( const BaseMeshPtr mesh ) const
{
    if ( !mesh )
        raiseAsterError( "ValueError: Mesh is empty" );

    if ( _mesh ){
        if( _mesh->getName() != mesh->getName() )
            raiseAsterError( "Incompatible meshes" );
    }
}

void Result::setMesh( const BaseMeshPtr &mesh )
{
    _checkMesh(mesh);
    _mesh = mesh;
};

void Result::setElementaryCharacteristics( const ElementaryCharacteristicsPtr &cara,
                                           ASTERINTEGER rank ) {

    if ( !cara )
        raiseAsterError( "ValueError: ElementaryCharacteristics is empty" );

    _mapElemCara[rank] = cara;
    ASTERINTEGER rang = rank;
    std::string type( "CARAELEM" );
    CALLO_RSADPA_ZK8_WRAP( getName(), &rang, cara->getName(), type );
    setMesh(cara->getMesh());
};

void Result::setListOfLoads( const ListOfLoadsPtr &load, ASTERINTEGER rank ) {
    _mapLoads[rank] = load;
    ASTERINTEGER rang = rank;
    std::string type( "EXCIT" );
    CALLO_RSADPA_ZK24_WRAP( getName(), &rang, load->getName(), type );
};

void Result::setMaterialField( const MaterialFieldPtr &mater, ASTERINTEGER rank ) {

    if ( !mater )
        raiseAsterError( "ValueError: MaterialField is empty" );

    _mapMaterial[rank] = mater;
    ASTERINTEGER rang = rank;
    std::string type( "CHAMPMAT" );
    CALLO_RSADPA_ZK8_WRAP( getName(), &rang, mater->getName(), type );
    setMesh( mater->getMesh() );
};

void Result::setModel( const ModelPtr &model, ASTERINTEGER rank ) {

    if ( !model )
        raiseAsterError( "ValueError: Model is empty" );

    _mapModel[rank] = model;
    ASTERINTEGER rang = rank;
    std::string type( "MODELE" );
    CALLO_RSADPA_ZK8_WRAP( getName(), &rang, model->getName(), type );
    const auto fed = model->getFiniteElementDescriptor();
    _fieldBuidler.addFiniteElementDescriptor( fed );
    setMesh(model->getMesh());
};

void Result::setTimeValue( ASTERDOUBLE value, ASTERINTEGER rank ) {
    ASTERINTEGER rang = rank;
    std::string type( "INST" );
    CALLO_RSADPA_ZR_WRAP( getName(), &rang, &value, type );
};

ASTERDOUBLE Result::getTimeValue( ASTERINTEGER rank ) {

    _serialNumber->updateValuePointer();
    _rspr->updateValuePointer();

    ASTERINTEGER nb_ranks = _serialNumber->size();

    AS_ASSERT( rank <= nb_ranks );

    auto &calcParam = _calculationParameter->getVectorOfObjects();
    auto nbParam = calcParam.size();

    for ( int i = 0; i < nbParam; ++i ) {
        const auto item = calcParam[i];
        auto typevar = trim( item[3].toString() );

        if ( typevar == "ACCES" ) {
            auto var_name = trim( _accessVariables->getStringFromIndex( i + 1 ) );
            if ( var_name == "INST" ) {
                auto nosuff = trim( item[0].toString() );
                auto ivar = std::stoi( trim( item[1].toString() ) );
                auto nmax = std::stoi( trim( item[2].toString() ) );

                AS_ASSERT( nosuff == ".RSPR" )

                auto index = nmax * rank + ivar - 1;

                return ( *_rspr )[index];
            }
        }
    }

    AS_ASSERT( false );
    return 0.0;
};

bool Result::allocate( ASTERINTEGER nbRanks ) {

    std::string base( JeveuxMemoryTypesNames[Permanent] );
    ASTERINTEGER nbordr = nbRanks;
    CALLO_RSCRSD( base, getName(), getType(), &nbordr );
    _nbRanks = nbRanks;

    AS_ASSERT( _calculationParameter->build( true ) );
    AS_ASSERT( _namesOfFields->build( true ) );

    return true;
};

void Result::setElementaryCharacteristics( const ElementaryCharacteristicsPtr &cara ) {
    _serialNumber->updateValuePointer();
    ASTERINTEGER nbRanks = _serialNumber->size();
    for ( ASTERINTEGER rank = 0; rank < nbRanks; ++rank ) {
        const ASTERINTEGER iordr = ( *_serialNumber )[rank];
        if ( _mapElemCara.find( iordr ) == _mapElemCara.end() )
            setElementaryCharacteristics( cara, iordr );
    }
};

void Result::setMaterialField( const MaterialFieldPtr &mater ) {
    _serialNumber->updateValuePointer();
    ASTERINTEGER nbRanks = _serialNumber->size();
    for ( ASTERINTEGER rank = 0; rank < nbRanks; ++rank ) {
        const ASTERINTEGER iordr = ( *_serialNumber )[rank];
        if ( _mapMaterial.find( iordr ) == _mapMaterial.end() )
            setMaterialField( mater, iordr );
    }
};

void Result::setModel( const ModelPtr &model ) {
    _serialNumber->updateValuePointer();
    ASTERINTEGER nbRanks = _serialNumber->size();
    for ( ASTERINTEGER rank = 0; rank < nbRanks; ++rank ) {
        const ASTERINTEGER iordr = ( *_serialNumber )[rank];
        if ( _mapModel.find( iordr ) == _mapModel.end() )
            setModel( model, iordr );
    }
};

BaseDOFNumberingPtr Result::getEmptyDOFNumbering() {
    std::string resuName( getName() );
    std::string name( "12345678.00000          " );
    ASTERINTEGER a = 10, b = 14;
    CALLO_GNOMSD( resuName, name, &a, &b );
    DOFNumberingPtr retour( new DOFNumbering( name.substr( 0, 14 ) ) );
    _listOfDOFNum.push_back( retour );
    return retour;
};

FieldOnNodesRealPtr Result::getEmptyFieldOnNodesReal( const std::string name,
                                                      const ASTERINTEGER rank ) {

    if ( rank > _nbRanks || rank <= 0 )
        raiseAsterError( "IndexError: Rank '" + std::to_string( rank ) + "' is out of range" );

    auto rschex = _getNewFieldName( name, rank );
    CALLO_RSNOCH( getName(), name, &rank );
    std::string bis( rschex.second.c_str(), 19 );
    FieldOnNodesRealPtr result = boost::make_shared< FieldOnNodesReal >( bis );

    auto curIter = _dictOfVectorOfFieldOnNodesReal.find( name );
    if ( curIter == _dictOfVectorOfFieldOnNodesReal.end() ) {
        _dictOfVectorOfFieldOnNodesReal[name] = VectorOfFieldOnNodesReal( _nbRanks );
    }
    _dictOfVectorOfFieldOnNodesReal[name][rank - 1] = result;
    return result;
};

#ifdef ASTER_HAVE_MPI
BaseDOFNumberingPtr Result::getEmptyParallelDOFNumbering() {
    std::string resuName( getName() );
    std::string name( "12345678.00000          " );
    ASTERINTEGER a = 10, b = 14;
    CALLO_GNOMSD( resuName, name, &a, &b );
    ParallelDOFNumberingPtr retour =
        boost::make_shared< ParallelDOFNumbering >( name.substr( 0, 14 ) );
    _listOfDOFNum.push_back( retour );
    return retour;
};
#endif /* ASTER_HAVE_MPI */

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
    auto curIter = _mapElemCara.find( rank );
    if ( curIter == _mapElemCara.end() )
        raiseAsterError( "IndexError: Rank '" + std::to_string( rank ) + "' not found" );
    return ( *curIter ).second;
};

bool Result::hasElementaryCharacteristics( ASTERINTEGER rank ) const {
    auto curIter = _mapElemCara.find( rank );
    if ( curIter == _mapElemCara.end() )
        return false;
    return true;
};

bool Result::hasElementaryCharacteristics( ) const {
    if ( _mapElemCara.empty() )
        return false;
    return true;
};

ListOfLoadsPtr Result::getListOfLoads( ASTERINTEGER rank ) const {
    auto curIter = _mapLoads.find( rank );
    if ( curIter == _mapLoads.end() )
        raiseAsterError( "IndexError: Rank '" + std::to_string( rank ) + "' not found" );
    return ( *curIter ).second;
};

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
    auto curIter = _mapMaterial.find( rank );
    if ( curIter == _mapMaterial.end() )
        raiseAsterError( "IndexError: Rank '" + std::to_string( rank ) + "' not found" );
    return ( *curIter ).second;
};

BaseMeshPtr Result::getMesh() const {
    if ( _mesh != nullptr )
        return _mesh;
    const auto model = getModel();
    if ( model != nullptr )
        return model->getMesh();
    return nullptr;
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

ModelPtr Result::getModel( ASTERINTEGER rank ) const {
    auto curIter = _mapModel.find( rank );
    if ( curIter == _mapModel.end() )
        raiseAsterError( "IndexError: Rank '" + std::to_string( rank ) + "' not found" );
    return ( *curIter ).second;
};

ASTERINTEGER Result::getNumberOfRanks() const { return _serialNumber->size(); };

VectorLong Result::getRanks() const {
    VectorLong v;
    _serialNumber->updateValuePointer();
    for ( ASTERINTEGER j = 0; j < _serialNumber->size(); ++j ) {
        v.push_back( ( *_serialNumber )[j] );
    }
    return v;
};

FieldOnCellsRealPtr Result::getFieldOnCellsReal( const std::string name,
                                                 const ASTERINTEGER rank ) const {
    if ( rank >= _nbRanks || rank < 0 )
        raiseAsterError( "IndexError: Rank '" + std::to_string( rank ) + "' is out of range" );

    auto curIter = _dictOfVectorOfFieldOnCellsReal.find( trim( name ) );
    if ( curIter == _dictOfVectorOfFieldOnCellsReal.end() )
        raiseAsterError( "ValueError: Field " + name + " unknown in the results container" );

    return curIter->second[rank];
};

ConstantFieldOnCellsChar16Ptr
Result::getConstantFieldOnCellsChar16( const std::string name, const ASTERINTEGER rank ) const {
    if ( rank >= _nbRanks || rank < 0 )
        raiseAsterError( "IndexError: Rank '" + std::to_string( rank ) + "' is out of range" );

    auto curIter = _dictOfVectorOfConstantFieldOnCellsChar16.find( trim( name ) );
    if ( curIter == _dictOfVectorOfConstantFieldOnCellsChar16.end() )
        raiseAsterError( "ValueError: Field " + name + " unknown in the results container" );
    return curIter->second[rank];
};

PyObject *Result::getAccessParameters() const {

    PyObject *returnDict = PyDict_New();
    std::string var_name, str_val, typevar, nosuff;
    int ivar, nmax, index;

    CALL_JEMARQ();
    AS_ASSERT( _serialNumber->updateValuePointer() );
    _rspr->updateValuePointer();
    _rspi->updateValuePointer();
    _rsp8->updateValuePointer();
    _rs16->updateValuePointer();
    _rs24->updateValuePointer();
    AS_ASSERT( _calculationParameter->build() );

    ASTERINTEGER nb_ranks = _serialNumber->size();

    var_name = "NUME_ORDRE";
    PyObject *listValues = PyList_New( nb_ranks );
    for ( int j = 0; j < nb_ranks; ++j ) {
        PyList_SetItem( listValues, j, PyLong_FromLong( ( *_serialNumber )[j] ) );
    }
    PyDict_SetItemString( returnDict, var_name.c_str(), listValues );
    Py_DECREF( listValues );

    for ( int i = 0; i < ( _calculationParameter->getVectorOfObjects() ).size(); ++i ) {
        const auto item = _calculationParameter->getVectorOfObjects()[i];
        typevar = trim( item[3].toString() );

        if ( typevar == "ACCES" ) {
            var_name = trim( _accessVariables->getStringFromIndex( i + 1 ) );
            nosuff = trim( item[0].toString() );
            ivar = std::stoi( trim( item[1].toString() ) );
            nmax = std::stoi( trim( item[2].toString() ) );

            PyObject *listValues = PyList_New( nb_ranks );

            if ( nosuff == ".RSPI" ) {
                for ( int j = 0; j < nb_ranks; ++j ) {
                    index = nmax * ( j ) + ivar - 1;
                    PyList_SetItem( listValues, j, PyLong_FromLong( ( *_rspi )[index] ) );
                }
            }

            else if ( nosuff == ".RSPR" ) {
                for ( int j = 0; j < nb_ranks; ++j ) {
                    index = nmax * ( j ) + ivar - 1;
                    PyList_SetItem( listValues, j, PyFloat_FromDouble( ( *_rspr )[index] ) );
                }
            }

            else {
                for ( int j = 0; j < nb_ranks; ++j ) {
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
                        PyList_SetItem( listValues, j, Py_None );
                    } else {
                        PyList_SetItem( listValues, j, PyUnicode_FromString( str_val.c_str() ) );
                    }
                }
            }
            PyDict_SetItemString( returnDict, var_name.c_str(), listValues );
            Py_DECREF( listValues );
        }
    }
    CALL_JEDEMA();
    return returnDict;
}

VectorString Result::getFieldsOnNodesNames() const {
    VectorString names;
    names.reserve( _dictOfVectorOfFieldOnNodesReal.size() );

    for ( auto &it : _dictOfVectorOfFieldOnNodesReal ) {
        std::string name = it.first;
        names.push_back( trim( name ) );
    }
    return names;
};

VectorString Result::getFieldsOnCellsNames() const {
    VectorString names;
    names.reserve( _dictOfVectorOfFieldOnCellsReal.size() );

    for ( auto &it : _dictOfVectorOfFieldOnCellsReal ) {
        std::string name = it.first;
        names.push_back( trim( name ) );
    }
    return names;
};

VectorString Result::getConstantFieldsOnCellsNames() const {
    VectorString names;
    names.reserve( _dictOfVectorOfConstantFieldOnCellsChar16.size() );

    for ( auto &it : _dictOfVectorOfConstantFieldOnCellsChar16 ) {
        std::string name = it.first;
        names.push_back( trim( name ) );
    }
    return names;
};

FieldOnNodesRealPtr Result::getFieldOnNodesReal( const std::string name,
                                                 const ASTERINTEGER rank ) const {

    if ( rank >= _nbRanks || rank < 0 )
        raiseAsterError( "IndexError: Rank '" + std::to_string( rank ) + "' is out of range" );

    auto curIter = _dictOfVectorOfFieldOnNodesReal.find( trim( name ) );
    if ( curIter == _dictOfVectorOfFieldOnNodesReal.end() )
        raiseAsterError( "ValueError: Field " + name + " unknown in the results container" );

    return curIter->second[rank];
};

bool Result::setField( const FieldOnNodesRealPtr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    CALL_JEMARQ();

    if ( !field )
        raiseAsterError( "ValueError: field is empty" );

    if ( rank >= _nbRanks || rank < 0 )
        raiseAsterError( "IndexError: Rank '" + std::to_string( rank ) + "' is out of range" );

    auto trim_name = trim( name );
    auto rschex = _getNewFieldName( trim_name, rank );
    AS_ASSERT( rschex.first <= 100 );

    CALLO_RSNOCH( getName(), name, &rank );
    std::string internalName( rschex.second.c_str(), 19 );
    FieldOnNodesRealPtr result = boost::make_shared< FieldOnNodesReal >( internalName, *field );

    auto curIter = _dictOfVectorOfFieldOnNodesReal.find( trim_name );
    if ( curIter == _dictOfVectorOfFieldOnNodesReal.end() ) {
        auto index = _symbolicNamesOfFields->getIndexFromString( trim_name );
        _dictOfVectorOfFieldOnNodesReal[trim_name] = VectorOfFieldOnNodesReal(
                            _nbRanks, FieldOnNodesRealPtr( nullptr ) );
    }else if(_dictOfVectorOfFieldOnNodesReal[trim_name].size() != _nbRanks){
         _dictOfVectorOfFieldOnNodesReal[trim_name].resize(
                            _nbRanks, FieldOnNodesRealPtr( nullptr ) );
    }

    _dictOfVectorOfFieldOnNodesReal[trim_name][rank] = result;

    CALL_JEDEMA();
    return true;
};

bool Result::setField( const FieldOnCellsRealPtr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    CALL_JEMARQ();

    if ( !field )
        raiseAsterError( "ValueError: field is empty" );

    if ( rank >= _nbRanks || rank < 0 )
        raiseAsterError( "IndexError: Rank '" + std::to_string( rank ) + "' is out of range" );

    auto trim_name = trim( name );
    auto rschex = _getNewFieldName( trim_name, rank );
    AS_ASSERT( rschex.first <= 100 );

    CALLO_RSNOCH( getName(), name, &rank );
    std::string internalName( rschex.second.c_str(), 19 );
    FieldOnCellsRealPtr result = boost::make_shared< FieldOnCellsReal >( internalName, *field );

    auto curIter = _dictOfVectorOfFieldOnCellsReal.find( trim_name );
    if ( curIter == _dictOfVectorOfFieldOnCellsReal.end() ) {
        auto index = _symbolicNamesOfFields->getIndexFromString( trim_name );
         _dictOfVectorOfFieldOnCellsReal[trim_name] = VectorOfFieldOnCellsReal(
                            _nbRanks, FieldOnCellsRealPtr( nullptr ) );
    }else if(_dictOfVectorOfFieldOnCellsReal[trim_name].size() != _nbRanks){
        _dictOfVectorOfFieldOnCellsReal[trim_name].resize(
                            _nbRanks, FieldOnCellsRealPtr( nullptr ) );
    }

    _dictOfVectorOfFieldOnCellsReal[trim_name][rank] = result;

    CALL_JEDEMA();
    return true;
};

bool Result::setField( const ConstantFieldOnCellsChar16Ptr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    CALL_JEMARQ();

    if ( !field )
        raiseAsterError( "ValueError: field is empty" );

    if ( rank >= _nbRanks || rank < 0 )
        raiseAsterError( "IndexError: Rank '" + std::to_string( rank ) + "' is out of range" );

    auto trim_name = trim( name );
    auto rschex = _getNewFieldName( trim_name, rank );
    AS_ASSERT( rschex.first <= 100 );

    CALLO_RSNOCH( getName(), name, &rank );
    std::string internalName( rschex.second.c_str(), 19 );
    auto result = boost::make_shared< ConstantFieldOnCellsChar16 >( internalName, *field );

    auto curIter = _dictOfVectorOfConstantFieldOnCellsChar16.find( trim_name );
    if ( curIter == _dictOfVectorOfConstantFieldOnCellsChar16.end() ) {
        auto index = _symbolicNamesOfFields->getIndexFromString( trim_name );
         _dictOfVectorOfConstantFieldOnCellsChar16[trim_name] = VectorOfConstantFieldOnCellsChar16(
                            _nbRanks, ConstantFieldOnCellsChar16Ptr( nullptr ) );
    }else if(_dictOfVectorOfConstantFieldOnCellsChar16[trim_name].size() != _nbRanks){
         _dictOfVectorOfConstantFieldOnCellsChar16[trim_name].resize(
                            _nbRanks, ConstantFieldOnCellsChar16Ptr( nullptr ) );
    }

    _dictOfVectorOfConstantFieldOnCellsChar16[trim_name][rank] = result;

    CALL_JEDEMA();
    return true;
};

void Result::printListOfFields() const {
    std::cout << "Content of DataStructure : ";
    for ( auto curIter : _dictOfVectorOfFieldOnNodesReal ) {
        std::cout << curIter.first << " - ";
    }
    for ( auto curIter : _dictOfVectorOfFieldOnCellsReal ) {
        std::cout << curIter.first << " - ";
    }
    for ( auto curIter : _dictOfVectorOfConstantFieldOnCellsChar16 ) {
        std::cout << curIter.first << " - ";
    }
    std::cout << std::endl;
}

void Result::printInfo() const {
    ASTERINTEGER umess( 6 );
    CALLO_RSINFO( getName(), &umess );
}

bool Result::printMedFile( const std::string fileName, std::string medName, bool local ) const {
    LogicalUnitFile a( fileName, Binary, New );
    ASTERINTEGER retour = a.getLogicalUnit();
    CommandSyntax cmdSt( "IMPR_RESU" );

    SyntaxMapContainer dict;
    dict.container["FORMAT"] = "MED";
    dict.container["UNITE"] = retour;

    if ( getMesh()->isParallel() )
    {
        dict.container["PROC0"] = "NON";
        if ( !local )
            dict.container["FICHIER_UNIQUE"] = "OUI";
    }
    else
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

    try {
        ASTERINTEGER op = 39;
        CALL_EXECOP( &op );
    } catch ( ... ) {
        AS_ASSERT( false );
    }

    return true;
};

bool Result::build() {
    CALL_JEMARQ();
    _serialNumber->updateValuePointer();

    AS_ASSERT( _calculationParameter->build( true ) );
    AS_ASSERT( _namesOfFields->build( true ) );

    const auto numberOfSerialNum = _serialNumber->size();
    _nbRanks = numberOfSerialNum;
    BaseMeshPtr curMesh( nullptr );
    const ASTERINTEGER iordr = ( *_serialNumber )[_nbRanks - 1];
    if ( _mapModel.find( iordr ) != _mapModel.end() ) {
        curMesh = _mapModel[iordr]->getMesh();
    } else if ( _mesh != nullptr ) {
        curMesh = _mesh;
    }

    ASTERINTEGER cmpt = 1;
    for ( const auto curIter : _namesOfFields->getVectorOfObjects() ) {
        auto nomSymb = trim( _symbolicNamesOfFields->getStringFromIndex( cmpt ) );
        AS_ASSERT( numberOfSerialNum <= curIter.size() );

        for ( ASTERINTEGER rank = 0; rank < numberOfSerialNum; ++rank ) {
            std::string name( trim( curIter[rank].toString() ) );
            if ( name != "" ) {
                CALL_JEMARQ();
                const std::string questi( "TYPE_CHAMP" );
                const std::string typeco( "CHAMP" );
                ASTERINTEGER repi = 0, ier = 0;
                JeveuxChar32 repk( " " );
                const std::string arret( "C" );

                CALLO_DISMOI( questi, name, typeco, &repi, repk, arret, &ier );
                const std::string resu( trim( repk.toString() ) );

                if ( resu == "NOEU" ) {
                    const auto &iterField = _dictOfVectorOfFieldOnNodesReal.find( nomSymb );
                    if ( iterField == _dictOfVectorOfFieldOnNodesReal.end() )
                        _dictOfVectorOfFieldOnNodesReal[nomSymb] = VectorOfFieldOnNodesReal(
                            numberOfSerialNum, FieldOnNodesRealPtr( nullptr ) );
                    else if ( ASTERINTEGER( iterField->second.size() ) != numberOfSerialNum ) {
                        iterField->second.resize( numberOfSerialNum,
                                                  FieldOnNodesRealPtr( nullptr ) );
                    }

                    ASTERINTEGER test2 = _dictOfVectorOfFieldOnNodesReal[nomSymb][rank].use_count();
                    if ( test2 == 0 ) {
                        AS_ASSERT( curMesh != nullptr );
                        FieldOnNodesRealPtr result =
                            _fieldBuidler.buildFieldOnNodes< ASTERDOUBLE >( name );
                        result->setMesh(curMesh);
                        _dictOfVectorOfFieldOnNodesReal[nomSymb][rank] = result;
                    }
                } else if ( resu == "ELEM" || resu == "ELNO" || resu == "ELGA" ) {
                    const auto &iterField = _dictOfVectorOfFieldOnCellsReal.find( nomSymb );
                    if ( iterField == _dictOfVectorOfFieldOnCellsReal.end() )
                        _dictOfVectorOfFieldOnCellsReal[nomSymb] = VectorOfFieldOnCellsReal(
                            numberOfSerialNum, FieldOnCellsRealPtr( nullptr ) );
                    else if ( ASTERINTEGER( iterField->second.size() ) != numberOfSerialNum ) {
                        iterField->second.resize( numberOfSerialNum,
                                                  FieldOnCellsRealPtr( nullptr ) );
                    }

                    ASTERINTEGER test2 = _dictOfVectorOfFieldOnCellsReal[nomSymb][rank].use_count();
                    if ( test2 == 0 ) {
                        AS_ASSERT( curMesh != nullptr );
                        FieldOnCellsRealPtr result =
                            _fieldBuidler.buildFieldOnCells< ASTERDOUBLE >( name, curMesh );
                        auto iterModel = _mapModel.find( ( *_serialNumber )[rank] );
                        if ( iterModel != _mapModel.end() )
                            if ( not( ( ( *iterModel ).second )->isEmpty() ) )
                                result->setModel( ( *iterModel ).second );
                            else if ( not( hasMultipleModel() ) ) {
                                ModelPtr curModel = getModel();
                                if ( not( curModel->isEmpty() ) )
                                    result->setModel( curModel );
                            }
                        _dictOfVectorOfFieldOnCellsReal[nomSymb][rank] = result;
                    }
                } else if ( resu == "CART" ) {
                    const auto &iterField =
                        _dictOfVectorOfConstantFieldOnCellsChar16.find( nomSymb );
                    if ( iterField == _dictOfVectorOfConstantFieldOnCellsChar16.end() )
                        _dictOfVectorOfConstantFieldOnCellsChar16[nomSymb] =
                            VectorOfConstantFieldOnCellsChar16(
                                numberOfSerialNum, ConstantFieldOnCellsChar16Ptr( nullptr ) );
                    else if ( ASTERINTEGER( iterField->second.size() ) != numberOfSerialNum ) {
                        iterField->second.resize( numberOfSerialNum,
                                                  ConstantFieldOnCellsChar16Ptr( nullptr ) );
                    }

                    ASTERINTEGER test2 =
                        _dictOfVectorOfConstantFieldOnCellsChar16[nomSymb][rank].use_count();
                    if ( test2 == 0 ) {
                        AS_ASSERT( curMesh != nullptr );
                        ConstantFieldOnCellsChar16Ptr result =
                            _fieldBuidler.buildConstantFieldOnCells< JeveuxChar16 >( name,
                                                                                     curMesh );
                        _dictOfVectorOfConstantFieldOnCellsChar16[nomSymb][rank] = result;
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

ASTERBOOL Result::resize(ASTERINTEGER nbRanks){
    try{
        _nbRanks = nbRanks;
        CALLO_RSAGSD( getName(), &nbRanks );
        return true;
    }catch(const std::exception& e){
        std::cout << e.what() << std::endl;
        return false;
    }

}
