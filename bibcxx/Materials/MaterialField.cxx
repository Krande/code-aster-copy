/**
 * @file MaterialField.cxx
 * @brief Implementation de MaterialField
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

#include "astercxx.h"
#include <stdexcept>
#include <typeinfo>

#include "Materials/MaterialField.h"
#include "Supervis/CommandSyntax.h"
#include "Utilities/SyntaxDictionary.h"

#include "Materials/MaterialFieldBuilder.h"

MaterialField::MaterialField( const std::string &name, const MeshPtr &mesh )
    : _mesh( mesh ), _model( nullptr ), DataStructure( name, 8, "CHAM_MATER" ),
      _listOfMaterials( ConstantFieldOnCellsChar8Ptr(
          new ConstantFieldOnCellsChar8( getName() + ".CHAMP_MAT ", mesh ) ) ),
      _listOfTemperatures( ConstantFieldOnCellsRealPtr(
          new ConstantFieldOnCellsReal( getName() + ".TEMPE_REF ", mesh ) ) ),
      _behaviourField( ConstantFieldOnCellsRealPtr(
          new ConstantFieldOnCellsReal( getName() + ".COMPOR ", mesh ) ) ),
      _cvrcNom( JeveuxVectorChar8( getName() + ".CVRCNOM" ) ),
      _cvrcGd( JeveuxVectorChar8( getName() + ".CVRCGD" ) ),
      _cvrcVarc( JeveuxVectorChar8( getName() + ".CVRCVARC" ) ),
      _cvrcCmp( JeveuxVectorChar8( getName() + ".CVRCCMP" ) ){};

MaterialField::MaterialField( const std::string &name, const SkeletonPtr &mesh )
    : _mesh( mesh ), _model( nullptr ), DataStructure( name, 8, "CHAM_MATER" ),
      _listOfMaterials( ConstantFieldOnCellsChar8Ptr(
          new ConstantFieldOnCellsChar8( getName() + ".CHAMP_MAT ", mesh ) ) ),
      _listOfTemperatures( ConstantFieldOnCellsRealPtr(
          new ConstantFieldOnCellsReal( getName() + ".TEMPE_REF ", mesh ) ) ),
      _behaviourField( ConstantFieldOnCellsRealPtr(
          new ConstantFieldOnCellsReal( getName() + ".COMPOR ", mesh ) ) ),
      _cvrcNom( JeveuxVectorChar8( getName() + ".CVRCNOM" ) ),
      _cvrcGd( JeveuxVectorChar8( getName() + ".CVRCGD" ) ),
      _cvrcVarc( JeveuxVectorChar8( getName() + ".CVRCVARC" ) ),
      _cvrcCmp( JeveuxVectorChar8( getName() + ".CVRCCMP" ) ){};

#ifdef ASTER_HAVE_MPI
MaterialField::MaterialField( const std::string &name, const ParallelMeshPtr &mesh )
    : _mesh( mesh ), _model( nullptr ), DataStructure( name, 8, "CHAM_MATER" ),
      _listOfMaterials( ConstantFieldOnCellsChar8Ptr(
          new ConstantFieldOnCellsChar8( getName() + ".CHAMP_MAT ", mesh ) ) ),
      _listOfTemperatures( ConstantFieldOnCellsRealPtr(
          new ConstantFieldOnCellsReal( getName() + ".TEMPE_REF ", mesh ) ) ),
      _behaviourField( ConstantFieldOnCellsRealPtr(
          new ConstantFieldOnCellsReal( getName() + ".COMPOR ", mesh ) ) ),
      _cvrcNom( JeveuxVectorChar8( getName() + ".CVRCNOM" ) ),
      _cvrcGd( JeveuxVectorChar8( getName() + ".CVRCGD" ) ),
      _cvrcVarc( JeveuxVectorChar8( getName() + ".CVRCVARC" ) ),
      _cvrcCmp( JeveuxVectorChar8( getName() + ".CVRCCMP" ) ){};
#endif /* ASTER_HAVE_MPI */

bool MaterialField::buildWithoutExternalStateVariables() {
    MaterialFieldBuilder::buildClass( *this );

    return true;
};

std::vector< MaterialPtr > MaterialField::getVectorOfMaterial() const {
    std::vector< MaterialPtr > toReturn;
    for ( const auto &curIter : _materialsFieldEntity )
        for ( const auto &curIter2 : curIter.first )
            toReturn.push_back( curIter2 );
    return toReturn;
};

std::vector< PartOfMaterialFieldPtr > MaterialField::getVectorOfPartOfMaterialField() const {
    std::vector< PartOfMaterialFieldPtr > toReturn;
    for ( const auto &curIter : _materialsFieldEntity ) {
        PartOfMaterialFieldPtr toPush(
            new PartOfMaterialField( curIter.first, curIter.second ) );
        toReturn.push_back( toPush );
    }
    return toReturn;
};

bool MaterialField::hasExternalStateVariables() const {
    return _cvrcVarc->exists();
};

bool MaterialField::hasExternalStateVariables( const std::string &name ) {
    if ( _cvrcVarc->exists() ) {
        _cvrcVarc->updateValuePointer();
        JeveuxChar8 toTest( name );
        auto size = _cvrcVarc->size();
        for ( int i = 0; i < size; ++i ) {
            if ( ( *_cvrcVarc )[i] == toTest )
                return true;
        }
    }

    return false;
};

void MaterialField::addBehaviourOnMesh( BehaviourDefinitionPtr &curBehav ) {
    _behaviours.push_back(
        listOfBehavAndGrpsValue( curBehav, MeshEntityPtr( new AllMeshEntities() ) ) );
}

void MaterialField::addBehaviourOnGroupOfCells( BehaviourDefinitionPtr &curBehav,
                                                     std::string nameOfGroup ) {
    if ( !_mesh )
        throw std::runtime_error( "Mesh is not defined" );
    if ( !_mesh->hasGroupOfCells( nameOfGroup ) )
        throw std::runtime_error( nameOfGroup + " not in mesh" );

    _behaviours.push_back(
        listOfBehavAndGrpsValue( curBehav, MeshEntityPtr( new GroupOfCells( nameOfGroup ) ) ) );
}

void MaterialField::addBehaviourOnCell( BehaviourDefinitionPtr &curBehav,
                                             std::string nameOfCell ) {
    if ( !_mesh )
        throw std::runtime_error( "Mesh is not defined" );

    _behaviours.push_back(
        listOfBehavAndGrpsValue( curBehav, MeshEntityPtr( new Cell( nameOfCell ) ) ) );
}

void MaterialField::addMaterialsOnMesh( std::vector< MaterialPtr > curMaters ) {
    _materialsFieldEntity.push_back(
        listOfMatsAndGrpsValue( curMaters, MeshEntityPtr( new AllMeshEntities() ) ) );
}

void MaterialField::addMaterialsOnMesh( MaterialPtr &curMater ) {
    addMaterialsOnMesh( (std::vector< MaterialPtr >){curMater} );
}

void MaterialField::addMaterialsOnGroupOfCells( std::vector< MaterialPtr > curMaters,
                                                     VectorString namesOfGroup ) {
    if ( !_mesh )
        throw std::runtime_error( "Mesh is not defined" );
    for ( const auto &nameOfGroup : namesOfGroup )
        if ( !_mesh->hasGroupOfCells( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );

    _materialsFieldEntity.push_back(
        listOfMatsAndGrpsValue( curMaters, MeshEntityPtr( new GroupOfCells( namesOfGroup ) ) ) );
}

void MaterialField::addMaterialsOnGroupOfCells( MaterialPtr &curMater,
                                                     VectorString namesOfGroup ) {
    addMaterialsOnGroupOfCells( (std::vector< MaterialPtr >){curMater}, namesOfGroup );
}

void MaterialField::addMaterialsOnCell( std::vector< MaterialPtr > curMaters,
                                             VectorString namesOfCells ) {
    if ( !_mesh )
        throw std::runtime_error( "Mesh is not defined" );

    _materialsFieldEntity.push_back(
        listOfMatsAndGrpsValue( curMaters, MeshEntityPtr( new Cell( namesOfCells ) ) ) );
}

void MaterialField::addMaterialsOnCell( MaterialPtr &curMater, VectorString namesOfCells ) {
    addMaterialsOnCell( (std::vector< MaterialPtr >){curMater}, namesOfCells );
}


void MaterialField::addExternalStateVariables( PyObject *keywords ) {

    // Check input PyObject
    if ( !PyDict_Check( keywords ) && !PyList_Check( keywords ) && !PyTuple_Check( keywords ) )
        throw std::runtime_error( "Unexpected value for 'AFFE_VARC'." );

    // Create syntax
    CommandSyntax cmdSt( "code_aster.Cata.Commons.c_affe_varc.C_AFFE_VARC_EXTE" );
    PyObject *kwfact = PyDict_New();
    PyDict_SetItemString( kwfact, "AFFE_VARC", keywords );
    cmdSt.define( kwfact );

    // Get objects to create
    std::string modelName;
    modelName = " ";
    if ( getModel() != NULL ) {
        modelName = getModel()->getName();
    }
    modelName.resize( 8, ' ' );
    std::string materialFieldName = getName();
    materialFieldName.resize( 8, ' ' );
    std::string meshName = getMesh()->getName();
    meshName.resize( 8, ' ' );

    // Add external state variables in material field
    CALLO_AFVARC(materialFieldName, meshName, modelName);
    CALLO_CMTREF(materialFieldName, meshName);

    Py_DECREF( kwfact );
};
