/**
 * @file LinearStaticAnalysis.cxx
 * @brief Fichier source contenant le source du solveur de mecanique statique
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

#include <stdexcept>
#include <chrono>

#include "Algorithms/GenericAlgorithm.h"
#include "Algorithms/StaticMechanicalAlgorithm.h"
#include "Numbering/DOFNumbering.h"
#include "Discretization/DiscreteProblem.h"
#include "Supervis/Exceptions.h"
#include "Analysis/LinearStaticAnalysis.h"
#include "Supervis/CommandSyntax.h"
#include "MemoryManager/deleteTemporaryObjects.h"

LinearStaticAnalysis::LinearStaticAnalysis(
    const ModelPtr &model, const MaterialFieldPtr &mater,
    const ElementaryCharacteristicsPtr &cara )
    : _model( model ), _materialField( mater ), _linearSolver( BaseLinearSolverPtr() ),
      _timeStep( boost::make_shared< TimeStepper >()  ), _sief_elga(true),
      _study( boost::make_shared< PhysicalProblem >( _model, _materialField, cara ) ) {
    _timeStep->setValues( VectorReal( 1, 0. ) );
};

void LinearStaticAnalysis::_computeStress( StaticMechanicalContext &ctx ) {
    auto start = std::chrono::high_resolution_clock::now();

    const auto &study = ctx.getDiscreteProblem()->getPhysicalProblem();
    const auto &model = study->getModel();
    const auto &mater = study->getMaterialField();
    const auto &load = study->getListOfLoads();
    const auto &cara = study->getElementaryCharacteristics();
    const auto &mateco = study->getCodedMaterial();

    bool l_sief_elga = _sief_elga;
    bool l_strx_elga = model->existsMultiFiberBeam();

    std::string caraName = "        ";
    if( cara )
        caraName = cara->getName();

    auto times = _timeStep->getValues();

    ASTERINTEGER nbrank = times.size();

    CALLO_COMPSTRESSFIELD(ctx.getResult()->getName(), model->getName(), mater->getName(),
                          mateco->getName(), caraName, load->getName(),
                          (ASTERLOGICAL *)&l_sief_elga, (ASTERLOGICAL *)&l_strx_elga,
                          &nbrank, times.data());

    auto finish = std::chrono::high_resolution_clock::now();
    ctx._timer["Post"] += std::chrono::duration<ASTERDOUBLE>(finish - start).count();
}

ElasticResultPtr LinearStaticAnalysis::execute( ElasticResultPtr resultC ) {

    _study->getCodedMaterial()->allocate(true);

    if ( !_timeStep )
        throw std::runtime_error( "No time list" );
    if ( _timeStep->size() == 0 )
        resultC->allocate( 1 );
    else
        resultC->allocate( _timeStep->size() );

    // Define the discrete problem
    DiscreteProblemPtr dProblem( boost::make_shared< DiscreteProblem >( _study ) );

    if ( _model->getMesh()->isParallel() ) {
        if ( !_linearSolver->isHPCCompliant() )
            throw std::runtime_error( "ParallelMesh not allowed with this linear solver" );
    }
    // Build the linear solver (sd_solver)
    _linearSolver->_commandName = "MECA_STATIQUE";
    if( _model->xfemPreconditioningEnable() ) _linearSolver->enableXfem();
    _linearSolver->build();

    BaseDOFNumberingPtr dofNum1;
#ifdef ASTER_HAVE_MPI
    if ( _model->getMesh()->isParallel() )
        dofNum1 = resultC->getEmptyParallelDOFNumbering();
    else
#endif /* ASTER_HAVE_MPI */
        dofNum1 = resultC->getEmptyDOFNumbering();
    dofNum1->setModel(_study->getModel());
    dofNum1->setListOfLoads(_study->getListOfLoads());
    _study->setDOFNumbering(dofNum1);
    _study->computeDOFNumbering( );

    StaticMechanicalContext currentContext( dProblem, _linearSolver, resultC );
    typedef Algorithm< TimeStepper, StaticMechanicalContext, StaticMechanicalAlgorithm >
        MSAlgo;
    MSAlgo::runAllStepsOverAlgorithm( *_timeStep, currentContext );

    // Destruct matrix
    const std::string matass_name = currentContext.getStiffnessMatrix()->getName();
    CALLO_DETMATRIX(matass_name);

    // Compute Post-process:
    _computeStress(currentContext);

    // Cleaning
    cleanJeveuxMemory();

    auto timer = currentContext.getTimer();

    std::cout << std::scientific
              << "Temps CPU consommé dans le calcul "
              <<   timer["Matrix"] + timer["Rhs"] + timer["Solve"] + timer["Post"] + timer["Facto"]
              << "s dont:" << std::endl;
    std::cout << std::scientific
              << "*Calcul et assemblage de la matrice en "
              <<   timer["Matrix"]  << "s" << std::endl;
    std::cout << std::scientific
              << "*Calcul et assemblage du second membre en "
              <<   timer["Rhs"]  << "s" << std::endl;
     std::cout << std::scientific
              << "*Factorisation de la matrice en "
              <<   timer["Facto"]  << "s" << std::endl;
    std::cout << std::scientific
              << "*Résolution du système linéaire en "
              <<   timer["Solve"]  << "s" << std::endl;
    std::cout << std::scientific
              << "*Post-traitements en "
              <<   timer["Post"]  << "s" << std::endl;
    return resultC;
};
