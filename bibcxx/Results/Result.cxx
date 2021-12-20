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

#include "aster_fort_ds.h"
#include "aster_fort_jeveux.h"
#include "aster_fort_superv.h"
#include "aster_fort_utils.h"
#include "aster_fort_calcul.h"

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

void
Result::addElementaryCharacteristics( const ElementaryCharacteristicsPtr &cara,
                                                        int rank ) {

    if( !cara )
      raiseAsterError( "ValueError: ElementaryCharacteristics is empty" );

    _mapElemCara[rank] = cara;
    ASTERINTEGER rang = rank;
    std::string type( "CARAELEM" );
    CALLO_RSADPA_ZK8_WRAP( getName(), &rang, cara->getName(), type );
    setMesh(cara->getMesh());
};

void Result::addListOfLoads( const ListOfLoadsPtr &load,
                                               int rank ) {
    _mapLoads[rank] = load;
    ASTERINTEGER rang = rank;
    std::string type( "EXCIT" );
    CALLO_RSADPA_ZK24_WRAP( getName(), &rang, load->getName(), type );
};

void Result::addMaterialField( const MaterialFieldPtr &mater,
                                                  int rank ) {

    if( !mater )
      raiseAsterError( "ValueError: MaterialField is empty" );

    _mapMaterial[rank] = mater;
    ASTERINTEGER rang = rank;
    std::string type( "CHAMPMAT" );
    CALLO_RSADPA_ZK8_WRAP( getName(), &rang, mater->getName(), type );
    setMesh( mater->getMesh() );
};

void Result::addModel( const ModelPtr &model,
                                         int rank ) {

    if( !model )
      raiseAsterError( "ValueError: Model is empty" );

    _mapModel[rank] = model;
    ASTERINTEGER rang = rank;
    std::string type( "MODELE" );
    CALLO_RSADPA_ZK8_WRAP( getName(), &rang, model->getName(), type );
    const auto fed = model->getFiniteElementDescriptor();
    _fieldBuidler.addFiniteElementDescriptor( fed );
    setMesh(model->getMesh());
};

void Result::addTimeValue( ASTERDOUBLE value, int rank ) {
    ASTERINTEGER rang = rank;
    std::string type( "INST" );
    CALLO_RSADPA_ZR_WRAP( getName(), &rang, &value, type );
};

bool Result::allocate( int nbRanks ) {
    std::string base( "G" );
    ASTERINTEGER nbordr = nbRanks;
    CALLO_RSCRSD( base, getName(), getType(), &nbordr );

    AS_ASSERT( _calculationParameter->build( true ) );
    AS_ASSERT( _namesOfFields->build( true ) );
    return true;
};

void Result::appendElementaryCharacteristicsOnAllRanks
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

void Result::appendMaterialFieldOnAllRanks( const MaterialFieldPtr &mater ) {
    _serialNumber->updateValuePointer();
    ASTERINTEGER nbRanks = _serialNumber->usedSize();
    for ( int rank = 0; rank < nbRanks; ++rank ) {
        const ASTERINTEGER iordr = ( *_serialNumber )[rank];
        if ( _mapMaterial.find( iordr ) == _mapMaterial.end() )
            addMaterialField( mater, iordr );
    }
};

void Result::appendModelOnAllRanks( const ModelPtr &model ) {
    _serialNumber->updateValuePointer();
    ASTERINTEGER nbRanks = _serialNumber->usedSize();
    for ( int rank = 0; rank < nbRanks; ++rank ) {
        const ASTERINTEGER iordr = ( *_serialNumber )[rank];
        if ( _mapModel.find( iordr ) == _mapModel.end() )
            addModel( model, iordr );
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

#ifdef ASTER_HAVE_MPI
BaseDOFNumberingPtr Result::getEmptyParallelDOFNumbering() {
    std::string resuName( getName() );
    std::string name( "12345678.00000          " );
    ASTERINTEGER a = 10, b = 14;
    CALLO_GNOMSD( resuName, name, &a, &b );
    ParallelDOFNumberingPtr retour( new ParallelDOFNumbering( name.substr( 0, 14 ) ) );
    _listOfDOFNum.push_back( retour );
    return retour;
};
#endif /* ASTER_HAVE_MPI */

std::vector< ElementaryCharacteristicsPtr >
Result::getAllElementaryCharacteristics() const
{
    return unique(_mapElemCara);
};

ElementaryCharacteristicsPtr Result::getElementaryCharacteristics() {
    const auto cara = getAllElementaryCharacteristics();
    AS_ASSERT(cara.size() <= 1);

    if( cara.size() == 1 )
        return cara[0];

    return ElementaryCharacteristicsPtr( nullptr );
};

ElementaryCharacteristicsPtr
Result::getElementaryCharacteristics( int rank ) {
    auto curIter = _mapElemCara.find( rank );
    if ( curIter == _mapElemCara.end() )
        throw std::runtime_error( "Rank not found" );
    return ( *curIter ).second;
};

bool Result::hasElementaryCharacteristics( ASTERINTEGER rank ) const {
    return _mapElemCara.count( rank ) > 0;
};

bool Result::hasElementaryCharacteristics() const { return !_mapElemCara.empty(); };

ListOfLoadsPtr Result::getListOfLoads( int rank ) {
    auto curIter = _mapLoads.find( rank );
    if ( curIter == _mapLoads.end() )
        throw std::runtime_error( "Rank not found" );
    return ( *curIter ).second;
};

std::vector< MaterialFieldPtr > Result::getMaterialFields() const
{
    return unique(_mapMaterial);
};

MaterialFieldPtr Result::getMaterialField() {
    const auto mate = getMaterialFields();
    AS_ASSERT(mate.size() <= 1);

    if( mate.size() == 1 )
        return mate[0];

    return MaterialFieldPtr( nullptr );
};

MaterialFieldPtr
Result::getMaterialField( int rank ) {
    auto curIter = _mapMaterial.find( rank );
    if ( curIter == _mapMaterial.end() )
        throw std::runtime_error( "Rank not found" );
    return ( *curIter ).second;
};

BaseMeshPtr Result::getMesh()
{
    if( _mesh != nullptr )
        return _mesh;
    const auto model = getModel();
    if( model != nullptr )
        return model->getMesh();
    return nullptr;
};

bool Result::hasMultipleModel()
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

bool Result::hasModel( const ASTERINTEGER &rank ) const { return _mapModel.count( rank ) > 0; }

bool Result::hasMaterialField( const ASTERINTEGER &rank ) const {
    return _mapMaterial.count( rank ) > 0;
}

std::vector< ModelPtr > Result::getModels() const
{
    return unique(_mapModel);
};

ModelPtr Result::getModel() {
    if ( hasMultipleModel() ){
        throw std::runtime_error( "Error: multiple models" );
    }

    const auto models = getModels();
    AS_ASSERT(models.size() <= 1 );

    if(models.size() == 1)
        return models[0];

    return ModelPtr( nullptr );

};

ModelPtr Result::getModel( int rank )
{
    auto curIter = _mapModel.find( rank );
    if ( curIter == _mapModel.end() )
        throw std::runtime_error( "Rank not found" );
    return ( *curIter ).second;
};

int Result::getNumberOfRanks() const
{
    return _serialNumber->usedSize();
};

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

PyObject *Result::getAccessParameters() const
{

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

    ASTERINTEGER nb_ranks = getNumberOfRanks();

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
                                                                     const ASTERINTEGER rank ) const
{

    return _dictOfMapOfFieldOnNodesReal.at( name ).at( rank );
};

FieldOnNodesComplexPtr Result::getFieldOnNodesComplex( const std::string name,
                                                       const ASTERINTEGER rank ) const {

    return _dictOfMapOfFieldOnNodesComplex.at( name ).at( rank );
};

bool Result::setField( const FieldOnNodesRealPtr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    CALL_JEMARQ();

    if ( !field )
        raiseAsterError( "ValueError: field is empty" );

    auto trim_name = trim( name );
    auto rschex = _getNewFieldName( trim_name, rank );
    AS_ASSERT( rschex.first <= 101 );

    if ( rschex.first == 101 )
        resize( std::max( (int)1, 2 * getNumberOfRanks() ) );

    CALLO_RSNOCH( getName(), name, &rank );
    std::string internalName( rschex.second.c_str(), 19 );
    FieldOnNodesRealPtr result = boost::make_shared< FieldOnNodesReal >( internalName, *field );

    _fieldBuidler.addFieldOnNodesDescription( result->getDescription() );

    if ( _dictOfMapOfFieldOnNodesReal.count( trim_name ) == 0 ) {
        _dictOfMapOfFieldOnNodesReal[trim_name] = MapOfFieldOnNodesReal();
    }

    _dictOfMapOfFieldOnNodesReal[trim_name][rank] = result;

    CALL_JEDEMA();
    return true;
};

bool Result::setField( const FieldOnNodesComplexPtr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    CALL_JEMARQ();

    if ( !field )
        raiseAsterError( "ValueError: field is empty" );

    auto trim_name = trim( name );
    auto rschex = _getNewFieldName( trim_name, rank );
    AS_ASSERT( rschex.first <= 101 );

    if ( rschex.first == 101 )
        resize( std::max( (int)1, 2 * getNumberOfRanks() ) );

    CALLO_RSNOCH( getName(), name, &rank );
    std::string internalName( rschex.second.c_str(), 19 );
    FieldOnNodesComplexPtr result =
        boost::make_shared< FieldOnNodesComplex >( internalName, *field );

    _fieldBuidler.addFieldOnNodesDescription( result->getDescription() );

    if ( _dictOfMapOfFieldOnNodesComplex.count( trim_name ) == 0 ) {
        _dictOfMapOfFieldOnNodesComplex[trim_name] = MapOfFieldOnNodesComplex();
    }

    _dictOfMapOfFieldOnNodesComplex[trim_name][rank] = result;

    CALL_JEDEMA();
    return true;
};

bool Result::setField( const FieldOnCellsRealPtr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    CALL_JEMARQ();

    if ( !field )
        raiseAsterError( "ValueError: field is empty" );

    auto trim_name = trim( name );
    auto rschex = _getNewFieldName( trim_name, rank );
    AS_ASSERT( rschex.first <= 101 );

    if ( rschex.first == 101 )
        resize( std::max( (int)1, 2 * getNumberOfRanks() ) );

    CALLO_RSNOCH( getName(), name, &rank );
    std::string internalName( rschex.second.c_str(), 19 );
    FieldOnCellsRealPtr result = boost::make_shared< FieldOnCellsReal >( internalName, *field );

    if ( _dictOfMapOfFieldOnCellsReal.count( trim_name ) == 0 ) {
        _dictOfMapOfFieldOnCellsReal[trim_name] = MapOfFieldOnCellsReal();
    }

    _dictOfMapOfFieldOnCellsReal[trim_name][rank] = result;

    CALL_JEDEMA();
    return true;
};

bool Result::setField( const FieldOnCellsComplexPtr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    CALL_JEMARQ();

    if ( !field )
        raiseAsterError( "ValueError: field is empty" );

    auto trim_name = trim( name );
    auto rschex = _getNewFieldName( trim_name, rank );
    AS_ASSERT( rschex.first <= 101 );

    if ( rschex.first == 101 )
        resize( std::max( (int)1, 2 * getNumberOfRanks() ) );

    CALLO_RSNOCH( getName(), name, &rank );
    std::string internalName( rschex.second.c_str(), 19 );
    FieldOnCellsComplexPtr result =
        boost::make_shared< FieldOnCellsComplex >( internalName, *field );

    if ( _dictOfMapOfFieldOnCellsComplex.count( trim_name ) == 0 ) {
        _dictOfMapOfFieldOnCellsComplex[trim_name] = MapOfFieldOnCellsComplex();
    }

    _dictOfMapOfFieldOnCellsComplex[trim_name][rank] = result;

    CALL_JEDEMA();
    return true;
};

bool Result::setField( const FieldOnCellsLongPtr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    CALL_JEMARQ();

    if ( !field )
        raiseAsterError( "ValueError: field is empty" );

    auto trim_name = trim( name );
    auto rschex = _getNewFieldName( trim_name, rank );
    AS_ASSERT( rschex.first <= 101 );

    if ( rschex.first == 101 )
        resize( std::max( (int)1, 2 * getNumberOfRanks() ) );

    CALLO_RSNOCH( getName(), name, &rank );
    std::string internalName( rschex.second.c_str(), 19 );
    FieldOnCellsLongPtr result = boost::make_shared< FieldOnCellsLong >( internalName, *field );

    if ( _dictOfMapOfFieldOnCellsLong.count( trim_name ) == 0 ) {
        _dictOfMapOfFieldOnCellsLong[trim_name] = MapOfFieldOnCellsLong();
    }

    _dictOfMapOfFieldOnCellsLong[trim_name][rank] = result;

    CALL_JEDEMA();
    return true;
};

bool Result::setField( const ConstantFieldOnCellsChar16Ptr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    CALL_JEMARQ();

    if ( !field )
        raiseAsterError( "ValueError: field is empty" );

    auto trim_name = trim( name );
    auto rschex = _getNewFieldName( trim_name, rank );
    AS_ASSERT( rschex.first <= 101 );

    if ( rschex.first == 101 )
        resize( std::max( (int)1, 2 * getNumberOfRanks() ) );

    CALLO_RSNOCH( getName(), name, &rank );
    std::string internalName( rschex.second.c_str(), 19 );
    auto result = boost::make_shared< ConstantFieldOnCellsChar16 >( internalName, *field );

    if ( _dictOfMapOfConstantFieldOnCellsChar16.count( trim_name ) == 0 ) {
        _dictOfMapOfConstantFieldOnCellsChar16[trim_name] = MapOfConstantFieldOnCellsChar16();
    }

    _dictOfMapOfConstantFieldOnCellsChar16[trim_name][rank] = result;

    CALL_JEDEMA();
    return true;
};

bool Result::setField( const ConstantFieldOnCellsRealPtr field, const std::string &name,
                       const ASTERINTEGER rank ) {
    CALL_JEMARQ();

    if ( !field )
        raiseAsterError( "ValueError: field is empty" );

    auto trim_name = trim( name );
    auto rschex = _getNewFieldName( trim_name, rank );
    AS_ASSERT( rschex.first <= 101 );

    if ( rschex.first == 101 )
        resize( std::max( (int)1, 2 * getNumberOfRanks() ) );

    CALLO_RSNOCH( getName(), name, &rank );
    std::string internalName( rschex.second.c_str(), 19 );
    auto result = boost::make_shared< ConstantFieldOnCellsReal >( internalName, *field );

    if ( _dictOfMapOfConstantFieldOnCellsReal.count( trim_name ) == 0 ) {
        _dictOfMapOfConstantFieldOnCellsReal[trim_name] = MapOfConstantFieldOnCellsReal();
    }

    _dictOfMapOfConstantFieldOnCellsReal[trim_name][rank] = result;

    CALL_JEDEMA();
    return true;
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

void Result::listFields() const {
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

bool Result::printMedFile( const std::string fileName,
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

bool Result::build()
{
    CALL_JEMARQ();
    _serialNumber->updateValuePointer();

    AS_ASSERT( _calculationParameter->build( true ) );
    AS_ASSERT( _namesOfFields->build( true ) );

    const auto nbRanks = getNumberOfRanks();

    ASTERINTEGER cmpt = 1;
    for ( const auto &curIter : _namesOfFields->getVectorOfObjects() ) {
        auto nomSymb = trim( _symbolicNamesOfFields->getStringFromIndex( cmpt ) );
        AS_ASSERT( nbRanks <= curIter.size() );

        for ( ASTERINTEGER index = 0; index < nbRanks; ++index ) {
            std::string name( trim( curIter[index].toString() ) );
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
                            FieldOnCellsRealPtr result =
                                _fieldBuidler.buildFieldOnCells< ASTERDOUBLE >( name, _mesh );
                            auto iterModel = _mapModel.find( rank );
                            if ( iterModel != _mapModel.end() ) {
                                if ( !( ( *iterModel ).second )->isEmpty() )
                                    result->setModel( ( *iterModel ).second );
                            } else if ( !hasMultipleModel() ) {
                                ModelPtr curModel = getModel();
                                if ( curModel && !curModel->isEmpty() )
                                    result->setModel( curModel );
                            }
                            _dictOfMapOfFieldOnCellsReal[nomSymb][rank] = result;
                        }
                    } else if ( scalaire == "C" ) {
                        if ( _dictOfMapOfFieldOnCellsComplex.count( nomSymb ) == 0 ) {
                            _dictOfMapOfFieldOnCellsComplex[nomSymb] = MapOfFieldOnCellsComplex();
                        }

                        if ( _dictOfMapOfFieldOnCellsComplex[nomSymb].count( rank ) == 0 ) {
                            AS_ASSERT( _mesh != nullptr );
                            FieldOnCellsComplexPtr result =
                                _fieldBuidler.buildFieldOnCells< ASTERCOMPLEX >( name, _mesh );
                            auto iterModel = _mapModel.find( rank );
                            if ( iterModel != _mapModel.end() ) {
                                if ( !( ( *iterModel ).second )->isEmpty() )
                                    result->setModel( ( *iterModel ).second );
                            } else if ( !hasMultipleModel() ) {
                                ModelPtr curModel = getModel();
                                if ( curModel && !curModel->isEmpty() )
                                    result->setModel( curModel );
                            }
                            _dictOfMapOfFieldOnCellsComplex[nomSymb][rank] = result;
                        }
                    } else if ( scalaire == "I" ) {
                        if ( _dictOfMapOfFieldOnCellsLong.count( nomSymb ) == 0 ) {
                            _dictOfMapOfFieldOnCellsLong[nomSymb] = MapOfFieldOnCellsLong();
                        }

                        if ( _dictOfMapOfFieldOnCellsLong[nomSymb].count( rank ) == 0 ) {
                            AS_ASSERT( _mesh != nullptr );
                            FieldOnCellsLongPtr result =
                                _fieldBuidler.buildFieldOnCells< ASTERINTEGER >( name, _mesh );
                            auto iterModel = _mapModel.find( rank );
                            if ( iterModel != _mapModel.end() ) {
                                if ( !( ( *iterModel ).second )->isEmpty() )
                                    result->setModel( ( *iterModel ).second );
                            } else if ( !hasMultipleModel() ) {
                                ModelPtr curModel = getModel();
                                if ( curModel && !curModel->isEmpty() )
                                    result->setModel( curModel );
                            }
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
                            ConstantFieldOnCellsChar16Ptr result =
                                _fieldBuidler.buildConstantFieldOnCells< JeveuxChar16 >( name,
                                                                                         _mesh );
                            _dictOfMapOfConstantFieldOnCellsChar16[nomSymb][rank] = result;
                        }
                    } else if ( scalaire == "R" ) {
                        if ( _dictOfMapOfConstantFieldOnCellsReal.count( nomSymb ) == 0 ) {
                            _dictOfMapOfConstantFieldOnCellsReal[nomSymb] =
                                MapOfConstantFieldOnCellsReal();
                        }

                        if ( _dictOfMapOfConstantFieldOnCellsReal[nomSymb].count( rank ) == 0 ) {
                            AS_ASSERT( _mesh != nullptr );
                            ConstantFieldOnCellsRealPtr result =
                                _fieldBuidler.buildConstantFieldOnCells< ASTERDOUBLE >( name,
                                                                                        _mesh );
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

bool Result::resize( ASTERINTEGER nbRanks ) {
    if ( getNumberOfRanks() == 0 ) {
        return allocate( nbRanks );
    } else {
        CALLO_RSAGSD( getName(), &nbRanks );
        return true;
    }

    return true;
}
