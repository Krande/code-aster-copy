/**
 * @file ResultBalancer.cxx
 * @brief Implementation de ResultBalancer
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

#include "Results/ResultBalancer.h"

#include "Materials/MaterialField.h"
#include "Meshes/MeshBalancer.h"
#include "ParallelUtilities/ArrayWrapper.h"
#include "ParallelUtilities/AsterMPI.h"

#include <set>

bool componentsCheck( const VectorString comp ) {
    VectorString compOut;
    AsterMPI::all_gather( comp, compOut );
    std::set< std::string > toCmp;
    for ( const auto &curCmp : comp )
        toCmp.insert( curCmp );
    for ( const auto &curCmp : compOut ) {
        if ( toCmp.count( curCmp ) == 0 ) {
            throw std::runtime_error(
                "Different components between processors is not yet allowed" );
        }
    }
    return true;
};

template < typename FieldType >
std::shared_ptr< SimpleFieldOnNodes< FieldType > >
balanceFieldOnNodesFromResult( const MeshBalancer &meshB, const ParallelMeshPtr newMesh,
                               const std::shared_ptr< SimpleFieldOnNodes< FieldType > > &field ) {
    field->updateValuePointers();
    typedef SimpleFieldOnNodes< FieldType > Field;
    const auto components = field->getComponents();
    componentsCheck( components );
    std::shared_ptr< Field > toReturn(
        new Field( newMesh, field->getPhysicalQuantity(), components, true ) );

    const auto &values = field->getValues();
    const auto &lValues = field->getLogicalValues();
    ( *toReturn )( 0, 0 ) = 3.0;
    return toReturn;
};

ResultPtr applyBalancingStrategy( const ResultPtr resuIn, const VectorInt &newLocalNodesList ) {
    auto meshIn = resuIn->getMesh();
    MeshBalancer meshB = MeshBalancer();
    meshB.buildFromBaseMesh( meshIn );
    auto meshOut = meshB.applyBalancingStrategy( newLocalNodesList );

    auto modelIn = resuIn->getModel();
    ModelPtr modelOut( new Model( meshOut, modelIn ) );
    modelOut->build();

    auto materIn = resuIn->getMaterialField();
    MaterialFieldPtr materOut( new MaterialField( meshOut, materIn ) );
    materOut->build();

    const auto &indexes = resuIn->getIndexes();
    const auto &lastIndex = resuIn->getLastIndex();

    const ListOfLoadsPtr lOLoadsIn = resuIn->getListOfLoads( lastIndex );
    ListOfLoadsPtr lOLoadsOut( new ListOfLoads( modelOut ) );
    const auto &BClist = lOLoadsIn->getDirichletBCs();
    for ( const auto &curBC : BClist ) {
        auto curBCOut = DirichletBCPtr( new DirichletBC( curBC, modelOut ) );
        curBCOut->buildFromSyntax();
        lOLoadsOut->addLoad( curBCOut );
    }

    ResultPtr resuOut( new Result( resuIn->getType() ) );
    resuOut->setModel( modelOut );
    resuOut->setMaterialField( materOut );
    resuOut->allocate( 1 );
    resuOut->setListOfLoads( lOLoadsIn, 1 );

    const auto fONList = resuIn->getFieldsOnNodesRealNames();
    for ( const auto &curName : fONList ) {
        std::cout << "Champ " << curName << std::endl;
        const auto &curFON = resuIn->getFieldOnNodesReal( curName, lastIndex );
        const auto &curSFON = toSimpleFieldOnNodes( curFON );
        auto newSFON = balanceFieldOnNodesFromResult( meshB, meshOut, curSFON );
        auto newFON = toFieldOnNodes( newSFON );
        resuOut->setField( newFON, curName, 1 );
    }

    return resuOut;
};
