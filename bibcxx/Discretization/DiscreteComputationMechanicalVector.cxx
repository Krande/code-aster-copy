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
                                            const FieldOnNodesRealPtr displ_step,
                                            const FieldOnCellsRealPtr stress,
                                            const FieldOnCellsRealPtr internVar,
                                            const ASTERDOUBLE &time_prev,
                                            const ASTERDOUBLE &time_step,
                                            const VectorString &groupOfCells ) const {

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
    CalculPtr calcul =
        createCalculForNonLinear( option, time_prev, time_prev + time_step, _externVarFieldPrev,
                                  _externVarFieldCurr, groupOfCells );
    FiniteElementDescriptorPtr FEDesc = calcul->getFiniteElementDescriptor();

    // Set current physical state
    calcul->addInputField( "PDEPLMR", displ );
    calcul->addInputField( "PDEPLPR", displ_step );
    calcul->addInputField( "PCONTMR", stress );
    calcul->addInputField( "PVARIMR", internVar );

    // Provisoire: pour TANGENTE=VERIFICATION, nécessité de variables internes à chaque itération
    FieldOnCellsRealPtr vari_iter =
        std::make_shared< FieldOnCellsReal >( FEDesc, currBehaviour, "ELGA_VARI_R", currElemChara );
    calcul->addInputField( "PVARIMP", vari_iter );

    // Create output vector
    auto elemVect = std::make_shared< ElementaryVectorReal >();
    elemVect->setModel( currModel );
    elemVect->setMaterialField( currMater );
    elemVect->setElementaryCharacteristics( currElemChara );
    elemVect->prepareCompute( option );

    // Create output fields
    FieldOnCellsRealPtr stress_curr =
        std::make_shared< FieldOnCellsReal >( FEDesc, nullptr, "ELGA_SIEF_R", currElemChara );
    FieldOnCellsLongPtr exitField = std::make_shared< FieldOnCellsLong >( FEDesc );
    FieldOnCellsRealPtr vari_curr =
        std::make_shared< FieldOnCellsReal >( FEDesc, currBehaviour, "ELGA_VARI_R", currElemChara );

    // Add output fields
    calcul->addOutputField( "PVARIPR", vari_curr );
    calcul->addOutputField( "PCONTPR", stress_curr );
    calcul->addOutputField( "PCODRET", exitField );

    // Add output elementary
    calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );

    // Compute
    FieldOnNodesRealPtr internalForces;
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) )
            elemVect->addElementaryTerm( calcul->getOutputElementaryTerm( "PVECTUR" ) );
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

/** @brief Compute AFFE_CHAR_MECA DDL_IMPO */
bool DiscreteComputation::addMecaImposedTerms( ElementaryVectorRealPtr elemVect,
                                               const ASTERDOUBLE time_value ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    // Main parameters
    const std::string calcul_option( "CHAR_MECA" );
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto listOfLoads = _phys_problem->getListOfLoads();
    elemVect->prepareCompute( calcul_option );

    // Get JEVEUX names of objects to call Fortran
    JeveuxVectorChar24 listOfLoadsList = listOfLoads->getListVector();
    JeveuxVectorLong listOfLoadsInfo = listOfLoads->getInformationVector();
    std::string nameLcha = ljust( listOfLoadsList->getName(), 24 );
    std::string nameInfc = ljust( listOfLoadsInfo->getName(), 24 );
    std::string vectElemName = ljust( elemVect->getName(), 24 );
    std::string modelName = ljust( currModel->getName(), 24 );
    std::string typres( "R" );

    CALLO_VEDIME( modelName, nameLcha, nameInfc, &time_value, typres, vectElemName );

    return true;
}

/** @brief Compute CHAR_MECA */
bool DiscreteComputation::addMecaNeumannTerms( ElementaryVectorRealPtr elemVect,
                                               const ASTERDOUBLE time_value,
                                               const ASTERDOUBLE time_delta,
                                               const ASTERDOUBLE time_theta,
                                               const FieldOnCellsRealPtr _externVarField ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    // Main parameters
    const std::string calcul_option( "CHAR_MECA" );
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto listOfLoads = _phys_problem->getListOfLoads();
    elemVect->prepareCompute( calcul_option );

    // Get JEVEUX names of objects to call Fortran
    JeveuxVectorChar24 listOfLoadsList = listOfLoads->getListVector();
    JeveuxVectorLong listOfLoadsInfo = listOfLoads->getInformationVector();
    std::string nameLcha = ljust( listOfLoadsList->getName(), 24 );
    std::string nameInfc = ljust( listOfLoadsInfo->getName(), 24 );
    std::string vectElemName = ljust( elemVect->getName(), 24 );
    std::string modelName = ljust( currModel->getName(), 24 );
    std::string materName = ljust( currMater->getName(), 24 );
    std::string currCodedMaterName =
        ljust( currCodedMater->getCodedMaterialField()->getName(), 24 );
    std::string stop( "S" );
    std::string currElemCharaName( " " );
    if ( currElemChara )
        currElemCharaName = currElemChara->getName();
    currElemCharaName.resize( 24, ' ' );

    // Get external state variables
    std::string externVarName( " " );
    if ( _externVarField ) {
        externVarName = _externVarField->getName();
    }
    externVarName.resize( 24, ' ' );
    CALLO_VECHME_WRAP( stop, modelName, nameLcha, nameInfc, &time_value, &time_delta, &time_theta,
                       currElemCharaName, materName, currCodedMaterName, vectElemName,
                       externVarName );

    return true;
}

FieldOnNodesRealPtr DiscreteComputation::contactForces( const MeshCoordinatesFieldPtr geom,
                                                        const FieldOnNodesRealPtr displ,
                                                        const FieldOnNodesRealPtr displ_step,
                                                        const FieldOnCellsRealPtr data ) const {
    // Select option for matrix
    std::string option = "CHAR_MECA_CONT";

    auto [Fed_Slave, Fed_pair] = _phys_problem->getListOfLoads()->getContactLoadDescriptor();

    // Prepare computing
    CalculPtr calcul = std::make_unique< Calcul >( option );
    calcul->setFiniteElementDescriptor( Fed_pair );

    // Set input field
    calcul->addInputField( "PGEOMER", _phys_problem->getMesh()->getCoordinates() );
    calcul->addInputField( "PGEOMCR", geom );
    calcul->addInputField( "PDEPL_M", displ );
    calcul->addInputField( "PDEPL_P", displ_step );
    calcul->addInputField( "PCONFR", data );

    // Create output vector
    auto elemVect = std::make_shared< ElementaryVectorDisplacementReal >();
    elemVect->prepareCompute( option );

    // Add output elementary
    calcul->addOutputElementaryTerm( "PVECTCR", std::make_shared< ElementaryTermReal >() );
    calcul->addOutputElementaryTerm( "PVECTFR", std::make_shared< ElementaryTermReal >() );

    // Computation
    calcul->compute();
    if ( calcul->hasOutputElementaryTerm( "PVECTCR" ) ) {
        elemVect->addElementaryTerm( calcul->getOutputElementaryTerm( "PVECTCR" ) );
    }
    if ( calcul->hasOutputElementaryTerm( "PVECTFR" ) ) {
        elemVect->addElementaryTerm( calcul->getOutputElementaryTerm( "PVECTFR" ) );
    }
    elemVect->build();

    return elemVect->assemble( _phys_problem->getDOFNumbering() );
}
