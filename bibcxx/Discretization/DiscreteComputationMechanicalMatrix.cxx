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

ElementaryMatrixDisplacementRealPtr DiscreteComputation::elasticStiffnessMatrix(
    const ASTERDOUBLE &time, const ASTERINTEGER &modeFourier, const VectorString &groupOfCells,
    const FieldOnCellsRealPtr _externVarField ) {

    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >();

    // Get main parameters
    const std::string option( "RIGI_MECA" );
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Check external state variables
    if ( currMater && currMater->hasExternalStateVariable() ) {
        if ( !_externVarField ) {
            AS_ABORT( "External state variables vector is missing" )
        }
    }

    // Set parameters of elementary matrix
    elemMatr->setModel( currModel );
    elemMatr->setMaterialField( currMater );
    elemMatr->setElementaryCharacteristics( currElemChara );

    // Check super-element
    if ( currModel->existsSuperElement() ) {
        std::string modelName = ljust( currModel->getName(), 8 );
        CALLO_CHECKSUPERELEMENT( option, modelName );
    }

    // Prepare computing
    CalculPtr _calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        _calcul->setModel( currModel );
    } else {
        _calcul->setGroupsOfCells( currModel, groupOfCells );
    }
    elemMatr->prepareCompute( option );

    // Add input fields
    _calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    if ( currMater ) {
        _calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
        _calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );
    }
    if ( _externVarField ) {
        _calcul->addInputField( "PVARCPR", _externVarField );
    }
    if ( currElemChara ) {
        _calcul->addElementaryCharacteristicsField( currElemChara );
    }
    _calcul->addFourierModeField( modeFourier );
    _calcul->addTimeField( "PTEMPSR", time );
    if ( currModel->existsXfem() ) {
        XfemModelPtr currXfemModel = currModel->getXfemModel();
        _calcul->addXFEMField( currXfemModel );
    }

    // Add output elementary terms
    _calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
    _calcul->addOutputElementaryTerm( "PMATUNS", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for rigidity
    if ( currModel->existsFiniteElement() ) {
        _calcul->compute();
        if ( _calcul->hasOutputElementaryTerm( "PMATUUR" ) )
            elemMatr->addElementaryTerm( _calcul->getOutputElementaryTerm( "PMATUUR" ) );
        if ( _calcul->hasOutputElementaryTerm( "PMATUNS" ) )
            elemMatr->addElementaryTerm( _calcul->getOutputElementaryTerm( "PMATUNS" ) );
    };

    // Compute elementary matrices for dual BC
    DiscreteComputation::baseDualStiffnessMatrix( _calcul, elemMatr );

    elemMatr->build();
    return elemMatr;
};
ElementaryMatrixDisplacementRealPtr
DiscreteComputation::massMatrix( const ASTERDOUBLE &time, const VectorString &groupOfCells,
                                 const FieldOnCellsRealPtr _externVarField ) {
    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >();
    // Get main parameters
    const std::string option( "MASS_MECA" );
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Check external state variables
    if ( currMater && currMater->hasExternalStateVariable() ) {
        if ( !_externVarField ) {
            AS_ABORT( "External state variables vector is missing" )
        }
    }

    // Set parameters of elementary matrix
    elemMatr->setModel( currModel );
    elemMatr->setMaterialField( currMater );
    elemMatr->setElementaryCharacteristics( currElemChara );

    // Check super-element
    if ( currModel->existsSuperElement() ) {
        std::string modelName = ljust( currModel->getName(), 8 );
        CALLO_CHECKSUPERELEMENT( option, modelName );
    }

    // Prepare computing
    auto _calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        _calcul->setModel( currModel );
    } else {
        _calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    elemMatr->prepareCompute( option );

    // Add input fields
    _calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    if ( currMater ) {
        _calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
        _calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );
    }
    if ( _externVarField ) {
        _calcul->addInputField( "PVARCPR", _externVarField );
    }
    if ( currElemChara ) {
        _calcul->addElementaryCharacteristicsField( currElemChara );
    }

    if ( currModel->existsXfem() ) {
        XfemModelPtr currXfemModel = currModel->getXfemModel();
        _calcul->addXFEMField( currXfemModel );
    }

    // Add output elementary terms
    _calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
    _calcul->addOutputElementaryTerm( "PMATUNS", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for mass
    if ( currModel->existsFiniteElement() ) {
        _calcul->compute();
        if ( _calcul->hasOutputElementaryTerm( "PMATUUR" ) )
            elemMatr->addElementaryTerm( _calcul->getOutputElementaryTerm( "PMATUUR" ) );
        if ( _calcul->hasOutputElementaryTerm( "PMATUNS" ) )
            elemMatr->addElementaryTerm( _calcul->getOutputElementaryTerm( "PMATUNS" ) );
    };

    elemMatr->build();
    return elemMatr;
};
ElementaryMatrixDisplacementRealPtr
DiscreteComputation::dampingMatrix( const ElementaryMatrixDisplacementRealPtr &massMatrix,
                                    const ElementaryMatrixDisplacementRealPtr &stiffnessMatrix,
                                    const ASTERDOUBLE &time, const VectorString &groupOfCells,
                                    const FieldOnCellsRealPtr _externVarField ) {
    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >();
    // Get main parameters
    const std::string option( "AMOR_MECA" );
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Check external state variables
    if ( currMater && currMater->hasExternalStateVariable() ) {
        if ( !_externVarField ) {
            AS_ABORT( "External state variables vector is missing" )
        }
    }

    // Set parameters of elementary matrix
    elemMatr->setModel( currModel );
    elemMatr->setMaterialField( currMater );
    elemMatr->setElementaryCharacteristics( currElemChara );

    // Check super-element
    if ( currModel->existsSuperElement() ) {
        std::string modelName = ljust( currModel->getName(), 8 );
        CALLO_CHECKSUPERELEMENT( option, modelName );
    }

    // Prepare computing
    auto _calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        _calcul->setModel( currModel );
    } else {
        _calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    elemMatr->prepareCompute( option );

    // Add input fields
    _calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    if ( currMater ) {
        _calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
        _calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );
    }
    if ( _externVarField ) {
        _calcul->addInputField( "PVARCPR", _externVarField );
    }
    if ( currElemChara ) {
        _calcul->addElementaryCharacteristicsField( currElemChara );
    }

    if ( massMatrix ) {
        ElementaryTermRealPtr resuElemMass;
        auto total_R_Mass = massMatrix->getElementaryTerms();
        for ( auto &R_Mass : total_R_Mass ) {
            if ( R_Mass->getFiniteElementDescriptor() == currModel->getFiniteElementDescriptor() ) {
                resuElemMass = R_Mass;
                break;
            }
        }
        _calcul->addInputElementaryTerm( "PMASSEL", resuElemMass );
    }

    if ( stiffnessMatrix ) {
        ElementaryTermRealPtr resuElemRigi;
        auto total_R_Rigi = stiffnessMatrix->getElementaryTerms();
        for ( auto &R_Rigi : total_R_Rigi ) {
            if ( R_Rigi->getFiniteElementDescriptor() == currModel->getFiniteElementDescriptor() ) {
                resuElemRigi = R_Rigi;
                break;
            }
        }

        if ( resuElemRigi->getPhysicalQuantity() == "MDNS_R" ) {
            _calcul->addInputElementaryTerm( "PRIGINS", resuElemRigi );
        } else {
            _calcul->addInputElementaryTerm( "PRIGIEL", resuElemRigi );
        }
    }

    // Add output elementary terms
    _calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
    _calcul->addOutputElementaryTerm( "PMATUNS", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for damping
    if ( currModel->existsFiniteElement() ) {
        _calcul->compute();
        if ( _calcul->hasOutputElementaryTerm( "PMATUUR" ) )
            elemMatr->addElementaryTerm( _calcul->getOutputElementaryTerm( "PMATUUR" ) );
        if ( _calcul->hasOutputElementaryTerm( "PMATUNS" ) )
            elemMatr->addElementaryTerm( _calcul->getOutputElementaryTerm( "PMATUNS" ) );
    };

    elemMatr->build();

    return elemMatr;
};
void DiscreteComputation::baseDualStiffnessMatrix( CalculPtr &calcul,
                                                   ElementaryMatrixDisplacementRealPtr &elemMatr ) {

    // Prepare loads
    const auto &_listOfLoads = _phys_problem->getListOfLoads();

    // Select option
    calcul->setOption( "MECA_DDLM_R" );
    calcul->clearInputs();
    calcul->clearOutputs();

    auto mecaLoadReal = _listOfLoads->getMechanicalLoadsReal();
    for ( const auto &curIter : mecaLoadReal ) {
        auto FEDesc = curIter->getFiniteElementDescriptor();
        auto field = curIter->getMechanicalLoadDescription()->getMultiplicativeField();
        if ( field && field->exists() && FEDesc ) {
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( FEDesc );
            calcul->addInputField( "PDDLMUR", field );
            calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATUUR" ) );
            }
        }
    }

    auto mecaLoadFunc = _listOfLoads->getMechanicalLoadsFunction();
    for ( const auto &curIter : mecaLoadFunc ) {
        auto FEDesc = curIter->getFiniteElementDescriptor();
        auto field = curIter->getMechanicalLoadDescription()->getMultiplicativeField();
        if ( field && field->exists() && FEDesc ) {
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( FEDesc );
            calcul->addInputField( "PDDLMUR", field );
            calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATUUR" ) );
            }
        }
    }

#ifdef ASTER_HAVE_MPI
    auto mecaParaLoadReal = _listOfLoads->getParallelMechanicalLoadsReal();
    for ( const auto &curIter : mecaParaLoadReal ) {
        auto FEDesc = curIter->getFiniteElementDescriptor();
        auto field = curIter->getMultiplicativeField();
        if ( field && field->exists() && FEDesc ) {
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( FEDesc );
            calcul->addInputField( "PDDLMUR", field );
            calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATUUR" ) );
            }
        }
    }

    auto mecaParaLoadFunc = _listOfLoads->getParallelMechanicalLoadsFunction();
    for ( const auto &curIter : mecaParaLoadFunc ) {
        auto FEDesc = curIter->getFiniteElementDescriptor();
        auto field = curIter->getMultiplicativeField();
        if ( field && field->exists() && FEDesc ) {
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( FEDesc );
            calcul->addInputField( "PDDLMUR", field );
            calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTerm( "PMATUUR" ) );
            }
        }
    }

#endif
};

ElementaryMatrixDisplacementRealPtr DiscreteComputation::dualStiffnessMatrix() {

    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >();

    // Get main parameters
    ModelPtr currModel = _phys_problem->getModel();
    MaterialFieldPtr currMater = _phys_problem->getMaterialField();
    ElementaryCharacteristicsPtr currElemChara = _phys_problem->getElementaryCharacteristics();

    // Set parameters of elementary matrix
    elemMatr->setModel( currModel );
    elemMatr->setMaterialField( currMater );
    elemMatr->setElementaryCharacteristics( currElemChara );

    // Prepare computing
    const std::string option( "MECA_DDLM_R" );
    CalculPtr _calcul = std::make_unique< Calcul >( option );
    elemMatr->prepareCompute( option );

    // Compute elementary matrices
    DiscreteComputation::baseDualStiffnessMatrix( _calcul, elemMatr );

    elemMatr->build();
    return elemMatr;
};

/** @brief Compute tangent matrix (not assembled) */
std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, ElementaryMatrixDisplacementRealPtr >
DiscreteComputation::computeTangentStiffnessMatrix(
    const FieldOnNodesRealPtr displ, const FieldOnNodesRealPtr displ_incr,
    const FieldOnCellsRealPtr stress, const FieldOnCellsRealPtr _internVar,
    const ConstantFieldOnCellsRealPtr _timeFieldPrev,
    const ConstantFieldOnCellsRealPtr _timeFieldCurr, const VectorString &groupOfCells ) {

    FieldOnCellsRealPtr _externVarFieldPrev;
    FieldOnCellsRealPtr _externVarFieldCurr;

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currBehaviour = _phys_problem->getBehaviourProperty();

    // Select option for matrix
    std::string option = "FULL_MECA";

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

    // Create output matrix
    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >();
    elemMatr->setModel( currModel );
    elemMatr->setMaterialField( currMater );
    elemMatr->setElementaryCharacteristics( currElemChara );
    elemMatr->prepareCompute( option );

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
    _calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
    _calcul->addOutputElementaryTerm( "PMATUNS", std::make_shared< ElementaryTermReal >() );
    _calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );

    // Compute
    if ( currModel->existsFiniteElement() ) {
        _calcul->compute();
        if ( _calcul->hasOutputElementaryTerm( "PMATUUR" ) )
            elemMatr->addElementaryTerm( _calcul->getOutputElementaryTerm( "PMATUUR" ) );
        if ( _calcul->hasOutputElementaryTerm( "PMATUNS" ) )
            elemMatr->addElementaryTerm( _calcul->getOutputElementaryTerm( "PMATUNS" ) );
        if ( _calcul->hasOutputElementaryTerm( "PVECTUR" ) )
            elemVect->addElementaryTerm( _calcul->getOutputElementaryTerm( "PVECTUR" ) );
        elemVect->build();
        elemMatr->build();
    };

    std::string exitFieldName = ljust( exitField->getName(), 19 );
    ASTERINTEGER exitCode = 0;
    CALLO_GETERRORCODE( exitFieldName, &exitCode );

#ifdef ASTER_HAVE_MPI
    ASTERINTEGER exitCodeLocal = exitCode;
    AsterMPI::all_reduce( exitCodeLocal, exitCode, MPI_MAX );
#endif

    return std::make_tuple( exitField, exitCode, elemMatr );
}

/** @brief Compute tangent prediction matrix (not assembled) */
std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, ElementaryMatrixDisplacementRealPtr >
DiscreteComputation::computeTangentPredictionMatrix(
    const FieldOnNodesRealPtr displ, const FieldOnNodesRealPtr displ_incr,
    const FieldOnCellsRealPtr stress, const FieldOnCellsRealPtr _internVar,
    const ConstantFieldOnCellsRealPtr _timeFieldPrev,
    const ConstantFieldOnCellsRealPtr _timeFieldCurr, const VectorString &groupOfCells ) {

    FieldOnCellsRealPtr _externVarFieldPrev;
    FieldOnCellsRealPtr _externVarFieldCurr;

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currBehaviour = _phys_problem->getBehaviourProperty();

    // Select option for matrix
    std::string option = "RIGI_MECA_TANG";

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

    // Create output matrix
    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >();
    elemMatr->setModel( currModel );
    elemMatr->setMaterialField( currMater );
    elemMatr->setElementaryCharacteristics( currElemChara );
    elemMatr->prepareCompute( option );

    // Create output vector
    auto elemVect = std::make_shared< ElementaryVectorDisplacementReal >();
    elemVect->setModel( currModel );
    elemVect->setMaterialField( currMater );
    elemVect->setElementaryCharacteristics( currElemChara );
    elemVect->prepareCompute( option );

    // Create output fields
    FieldOnCellsRealPtr stress_pred = std::make_shared< FieldOnCellsReal >(
        currModel, nullptr, "ELGA_SIEF_R", currElemChara, FEDesc = FEDesc );
    FieldOnCellsLongPtr maskField = std::make_shared< FieldOnCellsLong >( FEDesc );
    FieldOnCellsLongPtr exitField = std::make_shared< FieldOnCellsLong >( FEDesc );

    // Add output fields
    _calcul->addOutputField( "PCONTPR", stress_pred );
    _calcul->addOutputField( "PCOPRED", maskField );
    _calcul->addOutputField( "PCODRET", exitField );

    // Add output elementary
    _calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
    _calcul->addOutputElementaryTerm( "PMATUNS", std::make_shared< ElementaryTermReal >() );
    _calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );

    // Compute
    if ( currModel->existsFiniteElement() ) {
        _calcul->compute();
        if ( _calcul->hasOutputElementaryTerm( "PMATUUR" ) )
            elemMatr->addElementaryTerm( _calcul->getOutputElementaryTerm( "PMATUUR" ) );
        if ( _calcul->hasOutputElementaryTerm( "PMATUNS" ) )
            elemMatr->addElementaryTerm( _calcul->getOutputElementaryTerm( "PMATUNS" ) );
        if ( _calcul->hasOutputElementaryTerm( "PVECTUR" ) )
            elemVect->addElementaryTerm( _calcul->getOutputElementaryTerm( "PVECTUR" ) );
        elemVect->build();
        elemMatr->build();
    };

    std::string exitFieldName = ljust( exitField->getName(), 19 );
    ASTERINTEGER exitCode = 0;
    CALLO_GETERRORCODE( exitFieldName, &exitCode );
#ifdef ASTER_HAVE_MPI
    ASTERINTEGER exitCodeLocal = exitCode;
    AsterMPI::all_reduce( exitCodeLocal, exitCode, MPI_MAX );
#endif

    return std::make_tuple( exitField, exitCode, elemMatr );
}
