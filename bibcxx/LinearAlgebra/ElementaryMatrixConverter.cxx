/**
 * @file ElementaryMatrixConverter.cxx
 * @brief Implementation de ElementaryMatrixConverter
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

#include "LinearAlgebra/ElementaryMatrixConverter.h"

#ifdef ASTER_HAVE_MPI

ElementaryTermRealPtr transfertToParallelFEDesc( const ElementaryTermRealPtr curElemT,
                                                 const ParallelContactFEDescriptorPtr pFEDesc ) {
    const auto &fEDesc = pFEDesc->getSupportFiniteElementDescriptor();
    const auto &lielC = fEDesc->getListOfGroupsOfElementsExplorer();
    const auto &lielP = pFEDesc->getListOfGroupsOfElementsExplorer();
    const auto &cellMatching = pFEDesc->getCellMatching();

    const auto &values = curElemT->getValues();
    const auto name = values->getName();
    if ( !values->isBuilt() ) {
        values->build();
    }
    values->updateValuePointer();
    const auto &nbElemColl = values->size();

    ElementaryTermRealPtr termToReturn = std::make_shared< ElementaryTermReal >();
    termToReturn->allocate( pFEDesc, curElemT->getOption(), curElemT->getPhysicalQuantityId(),
                            curElemT->getLocalModeId() );

    auto &newVal = termToReturn->getValues();
    for ( int i = 0; i < nbElemColl; ++i ) {
        const auto &curLielP = lielP[i];
        const auto &cellNbC = lielC[i].getNumberOfNodes();
        const auto &cellNbP = curLielP.getNumberOfNodes();
        const auto &valueNb = ( *values )[i + 1]->size();
        const double matrixTermNbC = valueNb / cellNbC;
        const int valueNbP = matrixTermNbC * cellNbP;
        newVal->allocateObject( i + 1, valueNbP );
        ( *values )[i + 1]->updateValuePointer();
        ( *newVal )[i + 1]->updateValuePointer();
        int curPos = 0;
        for ( const auto &cellId : curLielP ) {
            const auto &cellIdC = cellMatching[i][curPos];
            for ( int iMatTerm = 0; iMatTerm < matrixTermNbC; ++iMatTerm ) {
                const auto pos1 = curPos * matrixTermNbC + iMatTerm;
                const auto pos2 = cellIdC * matrixTermNbC + iMatTerm;
                ( *( *newVal )[i + 1] )[pos1] = ( *( *values )[i + 1] )[pos2];
            }
            ++curPos;
        }
    }
    return termToReturn;
}

ElementaryMatrixDisplacementRealPtr
transfertToParallelFEDesc( const ElementaryMatrixDisplacementRealPtr eMIn,
                           const ParallelContactFEDescriptorPtr pFEDesc ) {
    auto eM = std::make_shared< ElementaryMatrixDisplacementReal >();
    auto elemTerms = eMIn->getElementaryTerms();
    const auto &fEDesc = pFEDesc->getSupportFiniteElementDescriptor();
    const auto &lielC = fEDesc->getListOfGroupsOfElementsExplorer();

#ifdef ASTER_DEBUG_CXX
    for ( const auto &colObj : lielC ) {
        for ( const auto &val : colObj ) {
            if ( val >= 0 )
                throw std::runtime_error( "Inconsistent situation" );
        }
    }
#endif

    for ( const auto &curElemT : elemTerms ) {
        const auto curFEDesc = curElemT->getFiniteElementDescriptor();
        if ( curFEDesc->getName() != fEDesc->getName() ) {
            throw std::runtime_error(
                "SupportFiniteElementDescriptor of "
                "ParallelContactFEDescriptor must be the same than  elementary matrix "
                "FiniteElementDescriptor" );
        }

        eM->addElementaryTerm( transfertToParallelFEDesc( curElemT, pFEDesc ) );
    }
    return eM;
};

ElementaryVectorDisplacementRealPtr
transfertToParallelFEDesc( const ElementaryVectorDisplacementRealPtr eVIn,
                           const ParallelContactFEDescriptorPtr pFEDesc ) {
    auto eV = std::make_shared< ElementaryVectorDisplacementReal >( eVIn->getModel() );
    auto elemTerms = eVIn->getElementaryTerms();
    const auto &fEDesc = pFEDesc->getSupportFiniteElementDescriptor();
    const auto &lielC = fEDesc->getListOfGroupsOfElementsExplorer();

#ifdef ASTER_DEBUG_CXX
    for ( const auto &colObj : lielC ) {
        for ( const auto &val : colObj ) {
            if ( val >= 0 )
                throw std::runtime_error( "Inconsistent situation" );
        }
    }
#endif

    for ( const auto &curElemT : elemTerms ) {
        const auto curFEDesc = curElemT->getFiniteElementDescriptor();
        if ( curFEDesc->getName() != fEDesc->getName() ) {
            throw std::runtime_error(
                "SupportFiniteElementDescriptor of "
                "ParallelContactFEDescriptor must be the same than  elementary matrix "
                "FiniteElementDescriptor" );
        }
        const auto newET = transfertToParallelFEDesc( curElemT, pFEDesc );
        newET->setFiniteElementDescriptor( pFEDesc );

        eV->addElementaryTerm( newET );
    }
    return eV;
};

#endif /* ASTER_HAVE_MPI */
