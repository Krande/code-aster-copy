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

#include "Results/Result.h"
#include "PythonBindings/LogicalUnitManager.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Tools.h"

void
ResultClass::addElementaryCharacteristics( const ElementaryCharacteristicsPtr &cara,
                                                        int rank ) {

    if( !cara )
      raiseAsterError( "ValueError: ElementaryCharacteristics is empty" );

    _mapElemCara[rank] = cara;
    ASTERINTEGER rang = rank;
    std::string type( "CARAELEM" );
    CALLO_RSADPA_ZK8_WRAP( getName(), &rang, cara->getName(), type );
};

void ResultClass::addListOfLoads( const ListOfLoadsPtr &load,
                                               int rank ) {
    _mapLoads[rank] = load;
    ASTERINTEGER rang = rank;
    std::string type( "EXCIT" );
    CALLO_RSADPA_ZK24_WRAP( getName(), &rang, load->getName(), type );
};

void ResultClass::addMaterialField( const MaterialFieldPtr &mater,
                                                  int rank ) {

    if( !mater )
      raiseAsterError( "ValueError: MaterialField is empty" );

    _mapMaterial[rank] = mater;
    ASTERINTEGER rang = rank;
    std::string type( "CHAMPMAT" );
    CALLO_RSADPA_ZK8_WRAP( getName(), &rang, mater->getName(), type );
};

void ResultClass::addModel( const ModelPtr &model,
                                         int rank ) {

    if( !model )
      raiseAsterError( "ValueError: Model is empty" );

    _mapModel[rank] = model;
    ASTERINTEGER rang = rank;
    std::string type( "MODELE" );
    CALLO_RSADPA_ZK8_WRAP( getName(), &rang, model->getName(), type );
    const auto fed = model->getFiniteElementDescriptor();
    _fieldBuidler.addFiniteElementDescriptor( fed );
};

void ResultClass::addTimeValue( ASTERDOUBLE value, int rank ) {
    ASTERINTEGER rang = rank;
    std::string type( "INST" );
    CALLO_RSADPA_ZR_WRAP( getName(), &rang, &value, type );
};

bool ResultClass::allocate( int nbRanks ) {
    std::string base( JeveuxMemoryTypesNames[getMemoryType()] );
    ASTERINTEGER nbordr = nbRanks;
    CALLO_RSCRSD( base, getName(), getType(), &nbordr );
    _nbRanks = nbRanks;
    return true;
};

void ResultClass::appendElementaryCharacteristicsOnAllRanks
    ( const ElementaryCharacteristicsPtr& cara )
{
    _serialNumber->updateValuePointer();
    ASTERINTEGER nbRanks = _serialNumber->usedSize();
    for ( int rank = 0; rank < nbRanks; ++rank ) {
        const ASTERINTEGER iordr = ( *_serialNumber )[rank];
        if ( _mapElemCara.find( iordr ) == _mapElemCara.end() )
            addElementaryCharacteristics( cara, iordr );
    }
};

void ResultClass::appendMaterialFieldOnAllRanks( const MaterialFieldPtr &mater ) {
    _serialNumber->updateValuePointer();
    ASTERINTEGER nbRanks = _serialNumber->usedSize();
    for ( int rank = 0; rank < nbRanks; ++rank ) {
        const ASTERINTEGER iordr = ( *_serialNumber )[rank];
        if ( _mapMaterial.find( iordr ) == _mapMaterial.end() )
            addMaterialField( mater, iordr );
    }
};

void ResultClass::appendModelOnAllRanks( const ModelPtr &model ) {
    _serialNumber->updateValuePointer();
    ASTERINTEGER nbRanks = _serialNumber->usedSize();
    for ( int rank = 0; rank < nbRanks; ++rank ) {
        const ASTERINTEGER iordr = ( *_serialNumber )[rank];
        if ( _mapModel.find( iordr ) == _mapModel.end() )
            addModel( model, iordr );
    }
};

BaseDOFNumberingPtr ResultClass::getEmptyDOFNumbering() {
    std::string resuName( getName() );
    std::string name( "12345678.00000          " );
    ASTERINTEGER a = 10, b = 14;
    CALLO_GNOMSD( resuName, name, &a, &b );
    DOFNumberingPtr retour( new DOFNumberingClass( name.substr( 0, 14 ) ) );
    _listOfDOFNum.push_back( retour );
    return retour;
};

FieldOnNodesRealPtr
ResultClass::getEmptyFieldOnNodesReal( const std::string name,
                                                      const int rank ) {

    if ( rank > _nbRanks || rank <= 0 )
      raiseAsterError( "IndexError: Rank '" + std::to_string(rank) + "' is out of range" );

    ASTERINTEGER retour;
    retour = 0;
    const ASTERINTEGER rankLong = rank;
    std::string null( " " );
    std::string returnName( 19, ' ' );
    CALLO_RSEXCH( null, getName(), name, &rankLong, returnName, &retour );
    CALLO_RSNOCH( getName(), name, &rankLong );
    std::string bis( returnName.c_str(), 19 );
    FieldOnNodesRealPtr result( new FieldOnNodesRealClass( bis ) );

    auto curIter = _dictOfVectorOfFieldsNodes.find( name );
    if ( curIter == _dictOfVectorOfFieldsNodes.end() ) {
        _dictOfVectorOfFieldsNodes[name] = VectorOfFieldsNodes( _nbRanks );
    }
    _dictOfVectorOfFieldsNodes[name][rank - 1] = result;
    return result;
};

#ifdef ASTER_HAVE_MPI
BaseDOFNumberingPtr ResultClass::getEmptyParallelDOFNumbering() {
    std::string resuName( getName() );
    std::string name( "12345678.00000          " );
    ASTERINTEGER a = 10, b = 14;
    CALLO_GNOMSD( resuName, name, &a, &b );
    ParallelDOFNumberingPtr retour( new ParallelDOFNumberingClass( name.substr( 0, 14 ) ) );
    _listOfDOFNum.push_back( retour );
    return retour;
};
#endif /* ASTER_HAVE_MPI */

std::vector< ElementaryCharacteristicsPtr >
ResultClass::getAllElementaryCharacteristics() const
{
    return unique(_mapElemCara);
};

ElementaryCharacteristicsPtr ResultClass::getElementaryCharacteristics() {
    const auto cara = getAllElementaryCharacteristics();
    AS_ASSERT(cara.size() <= 1);

    if( cara.size() == 1 )
        return cara[0];

    return ElementaryCharacteristicsPtr( nullptr );
};

ElementaryCharacteristicsPtr
ResultClass::getElementaryCharacteristics( int rank ) {
    auto curIter = _mapElemCara.find( rank );
    if ( curIter == _mapElemCara.end() )
        throw std::runtime_error( "Rank not found" );
    return ( *curIter ).second;
};

ListOfLoadsPtr ResultClass::getListOfLoads( int rank ) {
    auto curIter = _mapLoads.find( rank );
    if ( curIter == _mapLoads.end() )
        throw std::runtime_error( "Rank not found" );
    return ( *curIter ).second;
};

std::vector< MaterialFieldPtr > ResultClass::getMaterialFields() const
{
    return unique(_mapMaterial);
};

MaterialFieldPtr ResultClass::getMaterialField() {
    const auto mate = getMaterialFields();
    AS_ASSERT(mate.size() <= 1);

    if( mate.size() == 1 )
        return mate[0];

    return MaterialFieldPtr( nullptr );
};

MaterialFieldPtr
ResultClass::getMaterialField( int rank ) {
    auto curIter = _mapMaterial.find( rank );
    if ( curIter == _mapMaterial.end() )
        throw std::runtime_error( "Rank not found" );
    return ( *curIter ).second;
};

BaseMeshPtr ResultClass::getMesh()
{
    if( _mesh != nullptr )
        return _mesh;
    const auto model = getModel();
    if( model != nullptr )
        return model->getMesh();
    return nullptr;
};

bool ResultClass::hasMultipleModel()
{
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

std::vector< ModelPtr > ResultClass::getModels() const
{
    return unique(_mapModel);
};

ModelPtr ResultClass::getModel() {
    if ( hasMultipleModel() ){
        throw std::runtime_error( "Error: multiple models" );
    }

    const auto models = getModels();
    AS_ASSERT(models.size() <= 1 );

    if(models.size() == 1)
        return models[0];

    return ModelPtr( nullptr );

};

ModelPtr ResultClass::getModel( int rank )
{
    auto curIter = _mapModel.find( rank );
    if ( curIter == _mapModel.end() )
        throw std::runtime_error( "Rank not found" );
    return ( *curIter ).second;
};

int ResultClass::getNumberOfRanks() const
{
    return _serialNumber->usedSize();
};

VectorLong ResultClass::getRanks() const
{
    VectorLong v;
    _serialNumber->updateValuePointer();
    for ( int j = 0; j < _serialNumber->usedSize(); ++j ) {
        v.push_back( ( *_serialNumber )[j] );
    }
    return v;
};

FieldOnCellsRealPtr ResultClass::getFieldOnCellsReal( const std::string name,
                                                                           const int rank ) const
{
    if ( rank > _nbRanks || rank <= 0 )
      raiseAsterError( "IndexError: Rank '" + std::to_string(rank) + "' is out of range" );

    auto curIter = _dictOfVectorOfFieldsCells.find( trim( name ) );
    if ( curIter == _dictOfVectorOfFieldsCells.end() )
      raiseAsterError( "ValueError: Field " + name + " unknown in the results container" );

    FieldOnCellsRealPtr toReturn = curIter->second[rank - 1];
    return toReturn;
};

PyObject *ResultClass::getAccessParameters() const
{

  PyObject *returnDict = PyDict_New();
  std::string var_name, str_val, typevar, nosuff;
  int ivar, nmax, index;

  CALL_JEMARQ();
  _serialNumber->updateValuePointer();
  _rspr->updateValuePointer();
  _rspi->updateValuePointer();
  _rsp8->updateValuePointer();
  _rs16->updateValuePointer();
  _rs24->updateValuePointer();

  ASTERINTEGER nb_ranks = _serialNumber->usedSize();

  var_name = "NUME_ORDRE";
  PyObject *listValues = PyList_New( nb_ranks );
  for ( int j = 0; j < nb_ranks; ++j ) {
    PyList_SetItem(listValues, j, PyLong_FromLong( ( *_serialNumber )[j] ));
  }
  PyDict_SetItemString( returnDict, var_name.c_str(), listValues );
  Py_DECREF( listValues );

  for ( int i = 0; i < (_calculationParameter->getVectorOfObjects()).size(); ++i ) {
    const auto item = _calculationParameter->getVectorOfObjects()[i];
    typevar = trim(item[3].toString());

    if (typevar == "ACCES"){
      var_name = trim(_accessVariables->getStringFromIndex( i+1 ));
      nosuff = trim(item[0].toString());
      ivar = std::stoi(trim(item[1].toString()));
      nmax = std::stoi(trim(item[2].toString()));

      PyObject *listValues = PyList_New( nb_ranks );

      if (nosuff == ".RSPI") {
        for ( int j = 0; j < nb_ranks; ++j ) {
          index = nmax*(j)+ivar -1;
          PyList_SetItem( listValues, j, PyLong_FromLong( ( *_rspi )[index] ));
        }
      }

      else if (nosuff == ".RSPR") {
        for ( int j = 0; j < nb_ranks; ++j ) {
          index = nmax*(j)+ivar -1;
          PyList_SetItem( listValues, j, PyFloat_FromDouble( ( *_rspr )[index] ));
        }
      }

      else {
        for ( int j = 0; j < nb_ranks; ++j ) {
          index = nmax*(j)+ivar -1;
          if (nosuff == ".RSP8") {
            str_val = trim((( *_rsp8 )[index]).toString());
          }
          else if (nosuff == ".RS16") {
            str_val = trim((( *_rs16 )[index]).toString());
          }
          else if (nosuff == ".RS24") {
            str_val = trim((( *_rs24 )[index]).toString());
          }
          else {
            AS_ASSERT( false );
          }

          if (str_val.length()==0){
            PyList_SetItem( listValues, j, Py_None);
          }
          else {
            PyList_SetItem( listValues, j, PyUnicode_FromString( str_val.c_str()));
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

VectorString ResultClass::getFieldsOnNodesNames() const
{
  VectorString names;
  names.reserve( _dictOfVectorOfFieldsNodes.size());

  for ( auto& it : _dictOfVectorOfFieldsNodes ) {
    std::string name = it.first;
    names.push_back(trim(name)) ;
  }
  return names;
};

VectorString ResultClass::getFieldsOnCellsNames() const
{
  VectorString names;
  names.reserve( _dictOfVectorOfFieldsCells.size());

  for ( auto& it : _dictOfVectorOfFieldsCells ) {
    std::string name = it.first;
    names.push_back(trim(name)) ;
  }
  return names;
};


FieldOnNodesRealPtr ResultClass::getFieldOnNodesReal( const std::string name,
                                                                     const int rank ) const
{

    if ( rank > _nbRanks || rank <= 0 )
      raiseAsterError( "IndexError: Rank '" + std::to_string(rank) + "' is out of range" );

    auto curIter = _dictOfVectorOfFieldsNodes.find( trim( name ) );
    if ( curIter == _dictOfVectorOfFieldsNodes.end() )
      raiseAsterError( "ValueError: Field " + name + " unknown in the results container" );

    FieldOnNodesRealPtr toReturn = curIter->second[rank - 1];
    return toReturn;
};

void ResultClass::listFields() const
{
    std::cout << "Content of DataStructure : ";
    for ( auto curIter : _dictOfVectorOfFieldsNodes ) {
        std::cout << curIter.first << " - ";
    }
    for ( auto curIter : _dictOfVectorOfFieldsCells ) {
        std::cout << curIter.first << " - ";
    }
    std::cout << std::endl;
}

void ResultClass::printInfo() const {
    ASTERINTEGER umess( 6 );
    CALLO_RSINFO( getName(), &umess );
}

bool ResultClass::printMedFile( const std::string fileName,
                                std::string medName ) const
{
    LogicalUnitFile a( fileName, Binary, New );
    ASTERINTEGER retour = a.getLogicalUnit();
    CommandSyntax cmdSt( "IMPR_RESU" );

    SyntaxMapContainer dict;
    dict.container["FORMAT"] = "MED";
    dict.container["UNITE"] = retour;

    ListSyntaxMapContainer listeResu;
    SyntaxMapContainer dict2;
    dict2.container["RESULTAT"] = getName();
    dict2.container["TOUT_ORDRE"] = "OUI";
    if(!medName.empty())
        dict2.container["NOM_RESU_MED"] = medName.substr(0,8);
    listeResu.push_back( dict2 );
    dict.container["RESU"] = listeResu;

    cmdSt.define( dict );

    try {
        ASTERINTEGER op = 39;
        CALL_EXECOP( &op );
    } catch ( ... ) {
        AS_ASSERT(false);
    }

    return true;
};

bool ResultClass::update()
{
    CALL_JEMARQ();
    _serialNumber->updateValuePointer();

    AS_ASSERT( _calculationParameter->buildFromJeveux( true ) );
    AS_ASSERT( _namesOfFields->buildFromJeveux( true ) );

    const auto numberOfSerialNum = _serialNumber->usedSize();
    _nbRanks = numberOfSerialNum;
    BaseMeshPtr curMesh( nullptr );
    const ASTERINTEGER iordr = ( *_serialNumber )[_nbRanks - 1];
    if ( _mapModel.find( iordr ) != _mapModel.end() ){
        curMesh = _mapModel[iordr]->getMesh();
    }
    else if ( _mesh != nullptr ){
        curMesh = _mesh;
    }
    int cmpt = 1;
    for ( const auto curIter : _namesOfFields->getVectorOfObjects() ) {
        auto nomSymb = trim( _symbolicNamesOfFields->getStringFromIndex( cmpt ) );
        AS_ASSERT ( numberOfSerialNum <= curIter.size() );

        for ( int rank = 0; rank < numberOfSerialNum; ++rank ) {
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
                    const auto &iterField = _dictOfVectorOfFieldsNodes.find( nomSymb );
                    if ( iterField == _dictOfVectorOfFieldsNodes.end() )
                        _dictOfVectorOfFieldsNodes[nomSymb] = VectorOfFieldsNodes(
                            numberOfSerialNum, FieldOnNodesRealPtr( nullptr ) );
                    else if ( int(iterField->second.size()) != numberOfSerialNum ) {
                        iterField->second.resize( numberOfSerialNum,
                                                 FieldOnNodesRealPtr( nullptr ) );
                    }

                    ASTERINTEGER test2 = _dictOfVectorOfFieldsNodes[nomSymb][rank].use_count();
                    if ( test2 == 0 ) {
                        FieldOnNodesRealPtr result =
                            _fieldBuidler.buildFieldOnNodes< ASTERDOUBLE >( name );
                        _dictOfVectorOfFieldsNodes[nomSymb][rank] = result;
                    }
                } else if ( resu == "ELEM" || resu == "ELNO" || resu == "ELGA" ) {
                    const auto &iterField = _dictOfVectorOfFieldsCells.find( nomSymb );
                    if ( iterField == _dictOfVectorOfFieldsCells.end() )
                        _dictOfVectorOfFieldsCells[nomSymb] = VectorOfFieldsCells(
                            numberOfSerialNum, FieldOnCellsRealPtr( nullptr ) );
                    else if ( int(iterField->second.size()) != numberOfSerialNum ) {
                        iterField->second.resize( numberOfSerialNum,
                                                 FieldOnCellsRealPtr( nullptr ) );
                    }

                    ASTERINTEGER test2 = _dictOfVectorOfFieldsCells[nomSymb][rank].use_count();
                    if ( test2 == 0 ) {
                        AS_ASSERT( curMesh != nullptr );
                        FieldOnCellsRealPtr result =
                            _fieldBuidler.buildFieldOnCells< ASTERDOUBLE >( name, curMesh );
                        auto iterModel = _mapModel.find(( *_serialNumber )[rank]);
                        if ( iterModel != _mapModel.end() )
                            if ( not((( *iterModel ).second)->isEmpty()) )
                                result->setModel(( *iterModel ).second);
                        else if (not(hasMultipleModel())){
                            ModelPtr curModel = getModel();
                            if ( not(curModel->isEmpty()) )
                                result->setModel(curModel);
                        }
                        _dictOfVectorOfFieldsCells[nomSymb][rank] = result;

                    }
                }
                CALL_JEDEMA();
            }
        }
        ++cmpt;
    }

    CALL_JEDEMA();
    return update_tables();
};
