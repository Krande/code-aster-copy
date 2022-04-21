/**
 * @file MaterialField.cxx
 * @brief Implementation de MaterialField
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

#include "Materials/MaterialField.h"

#include "aster_fort_superv.h"
#include "astercxx.h"

#include "Modeling/FiniteElementDescriptor.h"
#include "Results/TransientResult.h"
#include "Supervis/CommandSyntax.h"
#include "Utilities/SyntaxDictionary.h"

MaterialField::MaterialField( const std::string &name, const MeshPtr &mesh )
    : _mesh( mesh ),
      _model( nullptr ),
      DataStructure( name, 8, "CHAM_MATER" ),
      _champ_mat( ConstantFieldOnCellsChar8Ptr(
          new ConstantFieldOnCellsChar8( getName() + ".CHAMP_MAT ", mesh ) ) ),
      _compor( ConstantFieldOnCellsRealPtr(
          new ConstantFieldOnCellsReal( getName() + ".COMPOR ", mesh ) ) ),
      _cvrcNom( JeveuxVectorChar8( getName() + ".CVRCNOM" ) ),
      _cvrcGd( JeveuxVectorChar8( getName() + ".CVRCGD" ) ),
      _cvrcVarc( JeveuxVectorChar8( getName() + ".CVRCVARC" ) ),
      _cvrcCmp( JeveuxVectorChar8( getName() + ".CVRCCMP" ) ){};

MaterialField::MaterialField( const std::string &name, const SkeletonPtr &mesh )
    : _mesh( mesh ),
      _model( nullptr ),
      DataStructure( name, 8, "CHAM_MATER" ),
      _champ_mat( ConstantFieldOnCellsChar8Ptr(
          new ConstantFieldOnCellsChar8( getName() + ".CHAMP_MAT ", mesh ) ) ),
      _compor( ConstantFieldOnCellsRealPtr(
          new ConstantFieldOnCellsReal( getName() + ".COMPOR ", mesh ) ) ),
      _cvrcNom( JeveuxVectorChar8( getName() + ".CVRCNOM" ) ),
      _cvrcGd( JeveuxVectorChar8( getName() + ".CVRCGD" ) ),
      _cvrcVarc( JeveuxVectorChar8( getName() + ".CVRCVARC" ) ),
      _cvrcCmp( JeveuxVectorChar8( getName() + ".CVRCCMP" ) ){};

#ifdef ASTER_HAVE_MPI
MaterialField::MaterialField( const std::string &name, const ParallelMeshPtr &mesh )
    : _mesh( mesh ),
      _model( nullptr ),
      DataStructure( name, 8, "CHAM_MATER" ),
      _champ_mat( ConstantFieldOnCellsChar8Ptr(
          new ConstantFieldOnCellsChar8( getName() + ".CHAMP_MAT ", mesh ) ) ),
      _compor( ConstantFieldOnCellsRealPtr(
          new ConstantFieldOnCellsReal( getName() + ".COMPOR ", mesh ) ) ),
      _cvrcNom( JeveuxVectorChar8( getName() + ".CVRCNOM" ) ),
      _cvrcGd( JeveuxVectorChar8( getName() + ".CVRCGD" ) ),
      _cvrcVarc( JeveuxVectorChar8( getName() + ".CVRCVARC" ) ),
      _cvrcCmp( JeveuxVectorChar8( getName() + ".CVRCCMP" ) ){};
#endif /* ASTER_HAVE_MPI */

listOfMaterials MaterialField::getVectorOfMaterial() const {
    listOfMaterials toReturn;
    for ( const auto &curIter : _materialsOnMeshEntities )
        for ( const auto &curIter2 : curIter.first )
            toReturn.push_back( curIter2 );
    return toReturn;
};

listOfPartOfMaterialField MaterialField::getVectorOfPartOfMaterialField() const {
    listOfPartOfMaterialField toReturn;
    for ( const auto &curIter : _materialsOnMeshEntities ) {
        PartOfMaterialFieldPtr toPush( new PartOfMaterialField( curIter.first, curIter.second ) );
        toReturn.push_back( toPush );
    }
    return toReturn;
};

void MaterialField::addBehaviourOnMesh( BehaviourDefinitionPtr &curBehav ) {
    _behaviourOnMeshEntities.push_back(
        listOfBehavioursOnMeshValue( curBehav, MeshEntityPtr( new AllMeshEntities() ) ) );
}

void MaterialField::addBehaviourOnGroupOfCells( BehaviourDefinitionPtr &curBehav,
                                                std::string nameOfGroup ) {
    if ( !_mesh )
        throw std::runtime_error( "Mesh is not defined" );
    if ( !_mesh->hasGroupOfCells( nameOfGroup ) )
        throw std::runtime_error( nameOfGroup + " not in mesh" );

    _behaviourOnMeshEntities.push_back(
        listOfBehavioursOnMeshValue( curBehav, MeshEntityPtr( new GroupOfCells( nameOfGroup ) ) ) );
}

void MaterialField::addMaterialsOnMesh( std::vector< MaterialPtr > curMaters ) {
    _materialsOnMeshEntities.push_back(
        listOfMaterialsOnMeshValue( curMaters, MeshEntityPtr( new AllMeshEntities() ) ) );
}

void MaterialField::addMaterialOnMesh( MaterialPtr &curMater ) {
    addMaterialsOnMesh( ( std::vector< MaterialPtr > ){ curMater } );
}

void MaterialField::addMaterialsOnGroupOfCells( std::vector< MaterialPtr > curMaters,
                                                VectorString namesOfGroup ) {
    if ( !_mesh )
        throw std::runtime_error( "Mesh is not defined" );
    for ( const auto &nameOfGroup : namesOfGroup )
        if ( !_mesh->hasGroupOfCells( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );

    _materialsOnMeshEntities.push_back( listOfMaterialsOnMeshValue(
        curMaters, MeshEntityPtr( new GroupOfCells( namesOfGroup ) ) ) );
}

void MaterialField::addMaterialOnGroupOfCells( MaterialPtr &curMater, VectorString namesOfGroup ) {
    addMaterialsOnGroupOfCells( ( std::vector< MaterialPtr > ){ curMater }, namesOfGroup );
}

ListSyntaxMapContainer MaterialField::syntaxForMaterial() {
    ListSyntaxMapContainer listAFFE;
    for ( auto &curIter : getMaterialsOnMeshEntities() ) {
        SyntaxMapContainer dict2;
        VectorString listOfMater;
        for ( const auto &curIter2 : curIter.first )
            listOfMater.push_back( curIter2->getName() );
        dict2.container["MATER"] = listOfMater;
        const MeshEntityPtr &tmp = curIter.second;
        if ( tmp->getType() == AllMeshEntitiesType )
            dict2.container["TOUT"] = "OUI";
        else if ( tmp->getType() == GroupOfCellsType )
            dict2.container["GROUP_MA"] = ( curIter.second )->getNames();
        else if ( tmp->getType() == GroupOfNodesType )
            dict2.container["GROUP_NO"] = ( curIter.second )->getNames();
        else
            throw std::runtime_error( "Support entity undefined" );
        listAFFE.push_back( dict2 );
    }
    return listAFFE;
}

ListSyntaxMapContainer MaterialField::syntaxForBehaviour() {
    ListSyntaxMapContainer listAFFE_COMPOR;
    for ( auto &curIter : getBehaviourOnMeshEntities() ) {
        SyntaxMapContainer dict2;
        dict2.container["COMPOR"] = curIter.first->getName();
        const MeshEntityPtr &tmp = curIter.second;
        if ( tmp->getType() == AllMeshEntitiesType )
            dict2.container["TOUT"] = "OUI";
        else if ( tmp->getType() == GroupOfCellsType )
            dict2.container["GROUP_MA"] = ( curIter.second )->getName();
        else
            throw std::runtime_error( "Support entity undefined" );
        listAFFE_COMPOR.push_back( dict2 );
    }
    return listAFFE_COMPOR;
}
