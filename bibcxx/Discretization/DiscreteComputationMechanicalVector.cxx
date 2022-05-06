/**
 * @file DiscreteComputation.cxx
 * @brief Implementation of class DiscreteComputation
 * @section LICENCE
 *   Copyright (C) 1991 2022  EDF R&D                www.code-aster.org
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

#include "aster_fort_calcul.h"
#include "aster_fort_superv.h"
#include "astercxx.h"

#include "Discretization/Calcul.h"
#include "Discretization/DiscreteComputation.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Modeling/XfemModel.h"
#include "Utilities/Tools.h"

/** @brief Compute internal forces, stress and internal state variables */
std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, FieldOnCellsRealPtr, FieldOnCellsRealPtr,
            FieldOnNodesRealPtr >
DiscreteComputation::computeInternalForces( const FieldOnNodesRealPtr displ,
                                            const FieldOnNodesRealPtr displ_incr,
                                            const FieldOnCellsRealPtr stress,
                                            const FieldOnCellsRealPtr _internVar,
                                            const ConstantFieldOnCellsRealPtr _timeFieldPrev,
                                            const ConstantFieldOnCellsRealPtr _timeFieldCurr,
                                            const VectorString &groupOfCells ) {

    FieldOnCellsRealPtr _externVarFieldPrev;
    FieldOnCellsRealPtr _externVarFieldCurr;

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currBehaviour = _phys_problem->getBehaviourProperty();

    // Select option for matrix
    std::string option = "RAPH_MECA";

    // Prepare computing:
    CalculPtr _calcul =
        createCalculForNonLinear( option, _timeFieldPrev, _timeFieldCurr, _externVarFieldPrev,
                                  _externVarFieldCurr, groupOfCells );
    FiniteElementDescriptorPtr FEDesc = _calcul->getFiniteElementDescriptor();

    // Set current physical state
    _calcul->addInputField( "PDEPLMR", displ );
    _calcul->addInputField( "PDEPLPR", displ_incr );
    _calcul->addInputField( "PCONTMR", stress );
    _calcul->addInputField( "PVARIMR", _internVar );

    // Provisoire: pour TANGENTE=VERIFICATION, nécessité de variables internes à chaque itération
    FieldOnCellsRealPtr vari_iter = std::make_shared< FieldOnCellsReal >(
        currModel, currBehaviour, "ELGA_VARI_R", currElemChara, FEDesc = FEDesc );
    _calcul->addInputField( "PVARIMP", vari_iter );

    // Create output vector
    auto elemVect = std::make_shared< ElementaryVectorDisplacementReal >();
    elemVect->setModel( currModel );
    elemVect->setMaterialField( currMater );
    elemVect->setElementaryCharacteristics( currElemChara );
    elemVect->prepareCompute( option );

    // Create output fields
    FieldOnCellsRealPtr stress_curr = std::make_shared< FieldOnCellsReal >(
        currModel, nullptr, "ELGA_SIEF_R", currElemChara, FEDesc = FEDesc );
    FieldOnCellsLongPtr exitField = std::make_shared< FieldOnCellsLong >( FEDesc );
    FieldOnCellsRealPtr vari_curr = std::make_shared< FieldOnCellsReal >(
        currModel, currBehaviour, "ELGA_VARI_R", currElemChara, FEDesc = FEDesc );

    // Add output fields
    _calcul->addOutputField( "PVARIPR", vari_curr );
    _calcul->addOutputField( "PCONTPR", stress_curr );
    _calcul->addOutputField( "PCODRET", exitField );

    // Add output elementary
    _calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );

    // Compute
    FieldOnNodesRealPtr internalForces;
    if ( currModel->existsFiniteElement() ) {
        _calcul->compute();
        if ( _calcul->hasOutputElementaryTerm( "PVECTUR" ) )
            elemVect->addElementaryTerm( _calcul->getOutputElementaryTerm( "PVECTUR" ) );
        elemVect->build();
        internalForces = elemVect->assemble( _phys_problem->getDOFNumbering() );
    };

    std::string exitFieldName = ljust( exitField->getName(), 19 );
    ASTERINTEGER exitCode = 0;
    CALLO_GETERRORCODE( exitFieldName, &exitCode );
#ifdef ASTER_HAVE_MPI
    ASTERINTEGER exitCodeLocal = exitCode;
    AsterMPI::all_reduce( exitCodeLocal, exitCode, MPI_MAX );
#endif

    return std::make_tuple( exitField, exitCode, vari_curr, stress_curr, internalForces );
}
