#ifndef RESULTBALANCER_H_
#define RESULTBALANCER_H_

/**
 * @file ResultBalancer.h
 * @brief Fichier entete de la classe ResultBalancer
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

#include "astercxx.h"

#include "Materials/MaterialField.h"
#include "Meshes/MeshBalancer.h"
#include "ParallelUtilities/ArrayWrapper.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Results/Result.h"

#ifdef ASTER_HAVE_MPI

bool componentsCheck( const VectorString comp );

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

    const auto compNb = components.size();
    auto values = field->getValues();
    auto lValues = field->getLogicalValues();
    auto wrapValues = ArrayWrapperPtr( new ArrayWrapper< JeveuxVectorReal >( values, compNb ) );
    JeveuxVectorReal medVB( "TemporaryObject" );
    ArrayWrapperPtr wrapB( new ArrayWrapper< JeveuxVectorReal >( medVB, compNb ) );
    auto nodeB = meshB.getNodeObjectBalancer();
    nodeB->balanceArrayOverProcessesWithRenumbering( wrapValues, wrapB );

    JeveuxVectorLong jvVL( "TemporaryObject3" );
    jvVL->allocate( lValues->size() );
    int cmpt = 0;
    for ( const auto &val : lValues ) {
        val ? ( *jvVL )[cmpt] = 1 : ( *jvVL )[cmpt] = 0;
        ++cmpt;
    }

    auto wrapLValues = ArrayWrapperLongPtr( new ArrayWrapper< JeveuxVectorLong >( jvVL, compNb ) );
    JeveuxVectorLong medVC( "TemporaryObject2" );
    ArrayWrapperLongPtr wrapC( new ArrayWrapper< JeveuxVectorLong >( medVC, compNb ) );
    nodeB->balanceArrayOverProcessesWithRenumbering( wrapLValues, wrapC );

    for ( int i = 0; i < newMesh->getNumberOfNodes(); ++i ) {
        for ( int j = 0; j < compNb; ++j ) {
            if ( ( *medVC )[i * compNb + j] == 1 ) {
                ( *toReturn )( i, j ) = ( *medVB )[i * compNb + j];
            }
        }
    }
    return toReturn;
};

template < typename FieldType >
std::shared_ptr< SimpleFieldOnCells< FieldType > >
balanceFieldOnCellsFromResult( const MeshBalancer &meshB, const ParallelMeshPtr newMesh,
                               const std::shared_ptr< SimpleFieldOnCells< FieldType > > &field ) {
    field->updateValuePointers();
    typedef SimpleFieldOnCells< FieldType > Field;
    const auto components = field->getComponents();
    componentsCheck( components );
    std::shared_ptr< Field > toReturn(
        new Field( newMesh, field->getLocalization(), field->getPhysicalQuantity(), components,
                   field->getMaxNumberOfPoints(), field->getMaxNumberOfSubPoints() ) );

    const auto compNb = components.size();
    const auto nbCells = field->getNumberOfCells();
    VectorLong cmps, lValues, cmps2;
    std::vector< FieldType > values;
    for ( int iCell = 0; iCell < nbCells; ++iCell ) {
        const auto nbPt = field->getNumberOfPointsOfCell( iCell );
        const auto nbSPt = field->getNumberOfSubPointsOfCell( iCell );
        const auto nbCmp = field->getNumberOfComponentsForSubpointsOfCell( iCell );
        for ( int iCmp = 0; iCmp < nbCmp; ++iCmp ) {
            for ( int iPt = 0; iPt < nbPt; ++iPt ) {
                for ( int iSPt = 0; iSPt < nbSPt; ++iSPt ) {
                    if ( field->hasValue( iCell, iCmp, iPt, iSPt ) ) {
                        values.push_back( field->getValue( iCell, iCmp, iPt, iSPt ) );
                        lValues.push_back( 1 );
                    } else {
                        values.push_back( 0. );
                        lValues.push_back( 0 );
                    }
                }
            }
        }
        cmps.push_back( nbPt * nbSPt * nbCmp );
        cmps2.push_back( nbPt );
        cmps2.push_back( nbSPt );
        cmps2.push_back( nbCmp );
    }

    auto wrapValues =
        ArrayWrapperVectorReal( new ArrayWrapper< std::vector< FieldType > >( values, cmps ) );
    VectorReal valuesOut;
    ArrayWrapperVectorReal wrapValueOut( new ArrayWrapper< VectorReal >( valuesOut, 0 ) );
    auto cellB = meshB.getCellObjectBalancer();
    cellB->balanceArrayOverProcessesWithRenumbering( wrapValues, wrapValueOut );

    auto wrapCmp2In = ArrayWrapperVectorLong( new ArrayWrapper< VectorLong >( cmps2, 3 ) );
    VectorLong cmps2Out;
    ArrayWrapperVectorLong wrapCmp2Out( new ArrayWrapper< VectorLong >( cmps2Out, 3 ) );
    cellB->balanceArrayOverProcessesWithRenumbering( wrapCmp2In, wrapCmp2Out );

    auto wrapLIn = ArrayWrapperVectorLong( new ArrayWrapper< VectorLong >( lValues, cmps ) );
    VectorLong lOut;
    ArrayWrapperVectorLong wrapLOut( new ArrayWrapper< VectorLong >( lOut, 0 ) );
    cellB->balanceArrayOverProcessesWithRenumbering( wrapLIn, wrapLOut );

    ASTERINTEGER cumSize = 0;
    const auto nbCellsNew = toReturn->getNumberOfCells();
    for ( int iCell = 0; iCell < nbCellsNew; ++iCell ) {
        const auto nbPt = cmps2Out[iCell * 3];
        const auto nbSPt = cmps2Out[iCell * 3 + 1];
        const auto nbCmp = cmps2Out[iCell * 3 + 2];
        for ( int iCmp = 0; iCmp < nbCmp; ++iCmp ) {
            for ( int iPt = 0; iPt < nbPt; ++iPt ) {
                for ( int iSPt = 0; iSPt < nbSPt; ++iSPt ) {
                    if ( lOut[cumSize] == 1 ) {
                        ( *toReturn )( iCell, iCmp, iPt, iSPt ) = valuesOut[cumSize];
                    }
                    ++cumSize;
                }
            }
        }
    }

    return toReturn;
};

/**
 * @fn applyBalancingStrategy
 * @brief Function to apply a balancing to a Result object
 * @author Nicolas Sellenet
 */
template < typename ResultType >
std::shared_ptr< ResultType > applyBalancingStrategy( const std::shared_ptr< ResultType > resuIn,
                                                      const VectorInt &newLocalNodesList ) {
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

    const auto &lastIndex = resuIn->getLastIndex();

    const ListOfLoadsPtr lOLoadsIn = resuIn->getListOfLoads( lastIndex );
    ListOfLoadsPtr lOLoadsOut( new ListOfLoads( modelOut ) );
    const auto &BClist = lOLoadsIn->getDirichletBCs();
    for ( const auto &curBC : BClist ) {
        if ( curBC->getType() == "CHAR_CINE_MECA" ) {
            auto curBC2 = std::static_pointer_cast< MechanicalDirichletBC >( curBC );
            auto curBCOut =
                MechanicalDirichletBCPtr( new MechanicalDirichletBC( curBC2, modelOut ) );
            curBCOut->buildFromSyntax();
            lOLoadsOut->addLoad( curBCOut );
        } else if ( curBC->getType() == "CHAR_CINE_THER" ) {
            auto curBC2 = std::static_pointer_cast< ThermalDirichletBC >( curBC );
            auto curBCOut = ThermalDirichletBCPtr( new ThermalDirichletBC( curBC2, modelOut ) );
            curBCOut->buildFromSyntax();
            lOLoadsOut->addLoad( curBCOut );
        } else if ( curBC->getType() == "CHAR_CINE_ACOU" ) {
            auto curBC2 = std::static_pointer_cast< AcousticDirichletBC >( curBC );
            auto curBCOut = AcousticDirichletBCPtr( new AcousticDirichletBC( curBC2, modelOut ) );
            curBCOut->buildFromSyntax();
            lOLoadsOut->addLoad( curBCOut );
        }
    }

    const auto &PMLlist = lOLoadsIn->getParallelMechanicalLoadsReal();
    for ( const auto &curBPML : PMLlist ) {
        auto curMLOut =
            ParallelMechanicalLoadRealPtr( new ParallelMechanicalLoadReal( curBPML, modelOut ) );
        lOLoadsOut->addLoad( curMLOut );
    }

    std::shared_ptr< ResultType > resuOut( new ResultType() );
    resuOut->allocate( 1 );

    const auto fONList = resuIn->getFieldsOnNodesRealNames();
    for ( const auto &curName : fONList ) {
        const auto &curFON = resuIn->getFieldOnNodesReal( curName, lastIndex );
        const auto &curSFON = toSimpleFieldOnNodes( curFON );
        auto newSFON = balanceFieldOnNodesFromResult( meshB, meshOut, curSFON );
        auto newFON = toFieldOnNodes( newSFON );
        resuOut->setField( newFON, curName, 1 );
    }

    const auto fOCList = resuIn->getFieldsOnCellsRealNames();
    const auto fED = modelOut->getFiniteElementDescriptor();
    for ( const auto &curName : fOCList ) {
        const auto &curFOC = resuIn->getFieldOnCellsReal( curName, lastIndex );
        const auto &curSFOC = toSimpleFieldOnCells( curFOC );
        auto newSFOC = balanceFieldOnCellsFromResult( meshB, meshOut, curSFOC );
        if ( meshIn->getNumberOfCells() != curSFOC->getNumberOfCells() )
            throw std::runtime_error( "Size inconsistency" );
        auto newFOC = toFieldOnCells( newSFOC, fED );
        resuOut->setField( newFOC, curName, 1 );
    }
    if ( resuIn->getFieldsOnNodesComplexNames().size() != 0 ||
         resuIn->getFieldsOnCellsComplexNames().size() != 0 ) {
        throw std::runtime_error( "Complex fields are not yet allowed" );
    }
    resuOut->setModel( modelOut );
    resuOut->setMaterialField( materOut );
    resuOut->setListOfLoads( lOLoadsOut, 1 );

    return resuOut;
};

#endif /* ASTER_HAVE_MPI */

#endif /* RESULTBALANCER_H_ */
