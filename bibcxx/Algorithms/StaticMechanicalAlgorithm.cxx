/**
 * @file StaticMechanicalAlgorithm.cxx
 * @brief Implementation de StaticMechanicalAlgorithm
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

#include "Algorithms/StaticMechanicalAlgorithm.h"

template <>
void updateContextFromStepper< TimeStepper::const_iterator, StaticMechanicalContext >(
    const TimeStepper::const_iterator &curStep, StaticMechanicalContext &context ) {
    context.setStep( *curStep, curStep.rank );
};

void StaticMechanicalAlgorithm::oneStep( const CurrentContext &ctx ) {
    BaseDOFNumberingPtr dofNum1 = ctx._results->getLastDOFNumbering();

    ctx._varCom->build( ctx._time );

    if ( ctx._rank == 1 || !ctx._isConst ) {
        auto matrElem = ctx._discreteProblem->computeElementaryStiffnessMatrix( ctx._time );

        // Build assembly matrix
        ctx._aMatrix->clearElementaryMatrix();
        ctx._aMatrix->appendElementaryMatrix( matrElem );
        ctx._aMatrix->setDOFNumbering( dofNum1 );
        ctx._aMatrix->setListOfLoads( ctx._listOfLoads );
        ctx._aMatrix->build();

        // Matrix factorization
        ctx._linearSolver->factorize( ctx._aMatrix );
    }

    // Build Dirichlet loads
    FieldOnNodesRealPtr chNoDir =
        ctx._discreteProblem->computeDirichlet( dofNum1, ctx._time );

    // Build Laplace forces
    ElementaryVectorPtr vectElem2 = ctx._discreteProblem->computeElementaryLaplaceVector();
    FieldOnNodesRealPtr chNoLap =
        vectElem2->assembleWithMultiplicatveFunction( dofNum1, ctx._time, Temporary );

    // Build Neumann loads
    VectorReal times;
    times.push_back( ctx._time );
    times.push_back( 0. );
    times.push_back( 0. );

    FieldOnNodesRealPtr chNoNeu =
        ctx._discreteProblem->computeNeumann( dofNum1, times, ctx._varCom );

    *chNoDir += *chNoLap;
    *chNoDir += *chNoNeu;

    if ( ctx._varCom->hasExternalStateVariables() ) {
        auto varComLoad = ctx._varCom->computeExternalStateVariablesLoad( dofNum1 );
        *chNoDir += *varComLoad;
    }

    CommandSyntax cmdSt( "MECA_STATIQUE" );
    cmdSt.setResult( ctx._results->getName(), ctx._results->getType() );

    FieldOnNodesRealPtr diriBCsFON =
        ctx._discreteProblem->computeDirichletBC( dofNum1, ctx._time, Temporary );

    FieldOnNodesRealPtr resultField =
        ctx._results->getEmptyFieldOnNodesReal( "DEPL", ctx._rank );

    resultField = ctx._linearSolver->solveWithDirichletBC(
        ctx._aMatrix, diriBCsFON, chNoDir, resultField );

    const auto &study = ctx._discreteProblem->getStudyDescription();
    const auto &model = study->getModel();
    const auto &mater = study->getMaterialField();
    const auto &load = study->getListOfLoads();
    const auto &cara = study->getElementaryCharacteristics();
    ctx._results->setModel( model, ctx._rank );
    ctx._results->setMaterialField( mater, ctx._rank );
    ctx._results->setTimeValue( ctx._time, ctx._rank );
    ctx._results->setListOfLoads( load, ctx._rank );
    if ( cara != nullptr ) {
        ctx._results->setElementaryCharacteristics( cara, ctx._rank );
    }
};
