/**
 * @file DiscreteComputation.cxx
 * @brief Implementation of class DiscreteComputation
 * @section LICENCE
 *   Copyright (C) 1991 2023  EDF R&D                www.code-aster.org
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

#include "astercxx.h"

#include "aster_fort_calcul.h"
#include "aster_fort_superv.h"

#include "Discretization/Calcul.h"
#include "Discretization/DiscreteComputation.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Messages/Messages.h"
#include "Modeling/Model.h"
#include "Modeling/XfemModel.h"
#include "Utilities/Tools.h"

ElementaryMatrixDisplacementRealPtr DiscreteComputation::getElasticStiffnessMatrix(
    const ASTERDOUBLE &time, const ASTERINTEGER &modeFourier, const VectorString &groupOfCells,
    const bool &with_dual ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    const std::string option( "RIGI_MECA" );

    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Check super-element
    if ( currModel->existsSuperElement() ) {
        std::string modelName = ljust( currModel->getName(), 8 );
        CALLO_CHECKSUPERELEMENT( option, modelName );
    }

    // Prepare computing
    CalculPtr calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
        calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );

        if ( currMater->hasExternalStateVariable() ) {
            calcul->addInputField( "PVARCPR", _phys_problem->getExternalStateVariables( time ) );
        }
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
    }
    calcul->addFourierModeField( modeFourier );
    calcul->addTimeField( "PTEMPSR", time );
    if ( currModel->existsXfem() ) {
        XfemModelPtr currXfemModel = currModel->getXfemModel();
        calcul->addXFEMField( currXfemModel );
    }

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
    calcul->addOutputElementaryTerm( "PMATUNS", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for rigidity
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUUR" ) );
        if ( calcul->hasOutputElementaryTerm( "PMATUNS" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUNS" ) );
    };

    if ( with_dual ) {
        DiscreteComputation::baseDualElasticStiffnessMatrix( calcul, elemMatr );
    }

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixDisplacementRealPtr DiscreteComputation::getGeometricStiffnessMatrix(
    const FieldOnCellsRealPtr sief_elga, const FieldOnCellsRealPtr strx_elga,
    const FieldOnNodesRealPtr displ, const ASTERINTEGER &modeFourier,
    const VectorString &groupOfCells ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    const std::string option( "RIGI_GEOM" );

    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Check super-element
    if ( currModel->existsSuperElement() ) {
        std::string modelName = ljust( currModel->getName(), 8 );
        CALLO_CHECKSUPERELEMENT( option, modelName );
    }

    // Prepare computing
    CalculPtr calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
    }

    calcul->addInputField( "PCONTRR", sief_elga );
    calcul->addInputField( "PEFFORR", sief_elga );

    if ( displ ) {
        calcul->addInputField( "PDEPLPR", displ );
    }

    if ( strx_elga ) {
        calcul->addInputField( "PSTRXRR", strx_elga );
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
    }

    calcul->addFourierModeField( modeFourier );

    if ( currModel->existsXfem() ) {
        XfemModelPtr currXfemModel = currModel->getXfemModel();
        calcul->addXFEMField( currXfemModel );
    }

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for rigidity
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUUR" ) );
    };

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixDisplacementRealPtr
DiscreteComputation::getFluidStructureStiffnessMatrix( const ASTERDOUBLE &time,
                                                       const ASTERINTEGER &modeFourier,
                                                       const VectorString &groupOfCells ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    const std::string option( "RIGI_FLUI_STRU" );

    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Check super-element
    if ( currModel->existsSuperElement() ) {
        std::string modelName = ljust( currModel->getName(), 8 );
        CALLO_CHECKSUPERELEMENT( option, modelName );
    }

    // Prepare computing
    CalculPtr calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
        calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );

        if ( currMater->hasExternalStateVariable() ) {
            calcul->addInputField( "PVARCPR", _phys_problem->getExternalStateVariables( time ) );
        }
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
    }
    calcul->addFourierModeField( modeFourier );

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for rigidity
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUUR" ) );
    };

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixDisplacementRealPtr
DiscreteComputation::getMechanicalMassMatrix( const bool diagonal, const ASTERDOUBLE &time,
                                              const VectorString &groupOfCells ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    std::string option( "MASS_MECA" );
    if ( diagonal )
        option = "MASS_MECA_DIAG";

    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Check super-element
    if ( currModel->existsSuperElement() ) {
        std::string modelName = ljust( currModel->getName(), 8 );
        CALLO_CHECKSUPERELEMENT( option, modelName );
    }

    // Prepare computing
    auto calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
        calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );

        if ( currMater->hasExternalStateVariable() ) {
            calcul->addInputField( "PVARCPR", _phys_problem->getExternalStateVariables( time ) );
        }
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
    }

    if ( currModel->existsXfem() ) {
        XfemModelPtr currXfemModel = currModel->getXfemModel();
        calcul->addXFEMField( currXfemModel );
    }

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
    if ( !diagonal ) {
        calcul->addOutputElementaryTerm( "PMATUNS", std::make_shared< ElementaryTermReal >() );
    }

    // Compute elementary matrices for mass
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUUR" ) );
        if ( !diagonal && calcul->hasOutputElementaryTerm( "PMATUNS" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUNS" ) );
    };

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixDisplacementRealPtr
DiscreteComputation::getFluidStructureMassMatrix( const ASTERDOUBLE &time,
                                                  const VectorString &groupOfCells ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    std::string option( "MASS_FLUI_STRU" );

    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Check super-element
    if ( currModel->existsSuperElement() ) {
        std::string modelName = ljust( currModel->getName(), 8 );
        CALLO_CHECKSUPERELEMENT( option, modelName );
    }

    // Prepare computing
    auto calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    calcul->addInputField( "PABSCUR", currModel->getMesh()->getCurvilinearAbscissa() );
    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
        calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );

        if ( currMater->hasExternalStateVariable() ) {
            calcul->addInputField( "PVARCPR", _phys_problem->getExternalStateVariables( time ) );
        }
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
    }

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for mass
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUUR" ) );
    };

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixDisplacementRealPtr DiscreteComputation::getMechanicalDampingMatrix(
    const ElementaryMatrixDisplacementRealPtr &massMatrix,
    const ElementaryMatrixDisplacementRealPtr &stiffnessMatrix, const ASTERDOUBLE &time,
    const VectorString &groupOfCells, const ASTERINTEGER &flui_int,
    const ASTERINTEGER &onde_flui ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    const std::string option( "AMOR_MECA" );
    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Check super-element
    if ( currModel->existsSuperElement() ) {
        std::string modelName = ljust( currModel->getName(), 8 );
        CALLO_CHECKSUPERELEMENT( option, modelName );
    }

    // Prepare computing
    auto calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    auto carteFluid = createDampingFluidField( flui_int, onde_flui );
    calcul->addInputField( "PAMORFL", carteFluid );
    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
        calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );

        if ( currMater->hasExternalStateVariable() ) {
            calcul->addInputField( "PVARCPR", _phys_problem->getExternalStateVariables( time ) );
        }
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
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
        calcul->addInputElementaryTerm( "PMASSEL", resuElemMass );
    }

    if ( stiffnessMatrix ) {
        std::vector< ElementaryTermRealPtr > resuElemRigi;
        auto total_R_Rigi = stiffnessMatrix->getElementaryTerms();
        for ( auto &R_Rigi : total_R_Rigi ) {
            if ( R_Rigi->getFiniteElementDescriptor() == currModel->getFiniteElementDescriptor() ) {
                resuElemRigi.push_back( R_Rigi );
            }
        }

        AS_ASSERT( resuElemRigi.size() <= 2 );
        if ( resuElemRigi[0]->getPhysicalQuantity() == "MDNS_R" ) {
            calcul->addInputElementaryTerm( "PRIGINS", resuElemRigi[0] );
            AS_ASSERT( resuElemRigi.size() == 1 );
        } else {
            calcul->addInputElementaryTerm( "PRIGIEL", resuElemRigi[0] );
            if ( resuElemRigi.size() > 1 ) {
                if ( resuElemRigi[1]->getPhysicalQuantity() == "MDNS_R" ) {
                    calcul->addInputElementaryTerm( "PRIGINS", resuElemRigi[1] );
                } else {
                    AS_ABORT( "Should be unsymetric" );
                }
            }
        }
    }

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
    calcul->addOutputElementaryTerm( "PMATUNS", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for damping
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUUR" ) );
        if ( calcul->hasOutputElementaryTerm( "PMATUNS" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUNS" ) );
    };

    elemMatr->build();

    return elemMatr;
};

ElementaryMatrixDisplacementComplexPtr DiscreteComputation::getHystereticStiffnessMatrix(
    const ElementaryMatrixDisplacementRealPtr &stiffnessMatrix, const ASTERDOUBLE &time,
    const VectorString &groupOfCells ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    const std::string option( "RIGI_MECA_HYST" );
    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementComplex >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Check super-element
    if ( currModel->existsSuperElement() ) {
        std::string modelName = ljust( currModel->getName(), 8 );
        CALLO_CHECKSUPERELEMENT( option, modelName );
    }

    // Prepare computing
    auto calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
        calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );

        if ( currMater->hasExternalStateVariable() ) {
            calcul->addInputField( "PVARCPR", _phys_problem->getExternalStateVariables( time ) );
        }
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
    }

    std::vector< ElementaryTermRealPtr > resuElemRigi;
    auto total_R_Rigi = stiffnessMatrix->getElementaryTerms();
    for ( auto &R_Rigi : total_R_Rigi ) {
        if ( R_Rigi->getFiniteElementDescriptor() == currModel->getFiniteElementDescriptor() ) {
            resuElemRigi.push_back( R_Rigi );
        }
    }

    AS_ASSERT( resuElemRigi.size() <= 2 );
    if ( resuElemRigi[0]->getPhysicalQuantity() == "MDNS_R" ) {
        calcul->addInputElementaryTerm( "PRIGINS", resuElemRigi[0] );
        AS_ASSERT( resuElemRigi.size() == 1 );
    } else {
        calcul->addInputElementaryTerm( "PRIGIEL", resuElemRigi[0] );
        if ( resuElemRigi.size() > 1 ) {
            if ( resuElemRigi[1]->getPhysicalQuantity() == "MDNS_R" ) {
                calcul->addInputElementaryTerm( "PRIGINS", resuElemRigi[1] );
            } else {
                AS_ABORT( "Should be unsymetric" );
            }
        }
    }

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PMATUUC", std::make_shared< ElementaryTermComplex >() );

    // Compute elementary matrices for complex rigidity
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATUUC" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermComplex( "PMATUUC" ) );
    };

    // This part is a hack because complex elem matr does not support real term
    // To remove later
    // Prepare loads
    const auto &listOfLoads = _phys_problem->getListOfLoads();

    // Select option
    calcul->setOption( "MECA_DDLM_RC" );

    auto impl = [&calcul, elemMatr]( auto loads ) {
        for ( const auto &load : loads ) {
            auto FEDesc = load->getFiniteElementDescriptor();
            auto field = load->getMultiplicativeField();
            if ( field && field->exists() && FEDesc && FEDesc->exists() ) {
                calcul->clearInputs();
                calcul->clearOutputs();
                calcul->setFiniteElementDescriptor( FEDesc );
                calcul->addInputField( "PDDLMUR", field );
                calcul->addOutputElementaryTerm( "PMATUUC",
                                                 std::make_shared< ElementaryTermComplex >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PMATUUC" ) ) {
                    elemMatr->addElementaryTerm(
                        calcul->getOutputElementaryTermComplex( "PMATUUC" ) );
                }
            }
        }
    };

    impl( listOfLoads->getMechanicalLoadsReal() );
    impl( listOfLoads->getMechanicalLoadsFunction() );
    impl( listOfLoads->getMechanicalLoadsComplex() );

#ifdef ASTER_HAVE_MPI
    impl( listOfLoads->getParallelMechanicalLoadsReal() );
    impl( listOfLoads->getParallelMechanicalLoadsFunction() );
#endif

    elemMatr->build();

    return elemMatr;
};

void DiscreteComputation::baseDualElasticStiffnessMatrix(
    CalculPtr &calcul, ElementaryMatrixDisplacementRealPtr &elemMatr ) const {

    // Prepare loads
    const auto &listOfLoads = _phys_problem->getListOfLoads();

    // Select option
    calcul->setOption( "MECA_DDLM_R" );

    auto impl = [calcul, elemMatr]( auto loads ) {
        for ( const auto &load : loads ) {
            auto FEDesc = load->getFiniteElementDescriptor();
            auto field = load->getMultiplicativeField();
            if ( field && field->exists() && FEDesc && FEDesc->exists() ) {
                calcul->clearInputs();
                calcul->clearOutputs();
                calcul->setFiniteElementDescriptor( FEDesc );
                calcul->addInputField( "PDDLMUR", field );
                calcul->addOutputElementaryTerm( "PMATUUR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) ) {
                    elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUUR" ) );
                }
            }
        }
    };

    // Real load
    impl( listOfLoads->getMechanicalLoadsReal() );
    impl( listOfLoads->getMechanicalLoadsFunction() );

#ifdef ASTER_HAVE_MPI
    impl( listOfLoads->getParallelMechanicalLoadsReal() );
    impl( listOfLoads->getParallelMechanicalLoadsFunction() );
#endif

    // Real part of complex load
    impl( listOfLoads->getMechanicalLoadsComplex() );
};

ElementaryMatrixDisplacementRealPtr DiscreteComputation::getDualElasticStiffnessMatrix() const {
    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    const std::string option( "MECA_DDLM_R" );

    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    // Prepare computing
    CalculPtr calcul = std::make_unique< Calcul >( option );

    // Compute elementary matrices
    DiscreteComputation::baseDualElasticStiffnessMatrix( calcul, elemMatr );

    elemMatr->build();
    return elemMatr;
};

/** @brief Compute tangent matrix (not assembled) */
std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, ElementaryMatrixDisplacementRealPtr >
DiscreteComputation::getTangentStiffnessMatrix(
    const FieldOnNodesRealPtr displ, const FieldOnNodesRealPtr displ_step,
    const FieldOnCellsRealPtr stress, const FieldOnCellsRealPtr internVar,
    const ASTERDOUBLE &time_prev, const ASTERDOUBLE &time_step,
    const FieldOnCellsRealPtr &externVarPrev, const FieldOnCellsRealPtr &externVarCurr,
    const VectorString &groupOfCells ) const {
    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currBehaviour = _phys_problem->getBehaviourProperty();

    // Select option for matrix
    std::string option = "FULL_MECA";

    // Prepare computing:
    CalculPtr calcul = createCalculForNonLinear( option, time_prev, time_prev + time_step,
                                                 externVarPrev, externVarCurr, groupOfCells );
    FiniteElementDescriptorPtr FEDesc = calcul->getFiniteElementDescriptor();

    // Set current physical state
    calcul->addInputField( "PDEPLMR", displ );
    calcul->addInputField( "PDEPLPR", displ_step );
    calcul->addInputField( "PCONTMR", stress );
    calcul->addInputField( "PVARIMR", internVar );

    // Provisoire: pour TANGENTE=VERIFICATION, nécessité de variables internes à chaque itération
    FieldOnCellsRealPtr vari_iter = std::make_shared< FieldOnCellsReal >(
        FEDesc, "ELGA", "VARI_R", currBehaviour, currElemChara );
    calcul->addInputField( "PVARIMP", internVar );

    // Create output matrix
    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    // Create output vector
    auto elemVect = std::make_shared< ElementaryVectorReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );
    elemVect->prepareCompute( option );

    // Create output fields
    FieldOnCellsRealPtr stress_curr =
        std::make_shared< FieldOnCellsReal >( FEDesc, "ELGA", "SIEF_R", currElemChara );
    FieldOnCellsLongPtr exitField = std::make_shared< FieldOnCellsLong >( FEDesc );
    FieldOnCellsRealPtr vari_curr = std::make_shared< FieldOnCellsReal >(
        FEDesc, "ELGA", "VARI_R", currBehaviour, currElemChara );

    // Add output fields
    calcul->addOutputField( "PVARIPR", vari_curr );
    calcul->addOutputField( "PCONTPR", stress_curr );
    calcul->addOutputField( "PCODRET", exitField );

    // Add output elementary
    calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
    calcul->addOutputElementaryTerm( "PMATUNS", std::make_shared< ElementaryTermReal >() );
    calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );

    // Compute
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUUR" ) );
        if ( calcul->hasOutputElementaryTerm( "PMATUNS" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUNS" ) );
        if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) )
            elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ) );
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
DiscreteComputation::getPredictionTangentStiffnessMatrix(
    const FieldOnNodesRealPtr displ, const FieldOnNodesRealPtr displ_step,
    const FieldOnCellsRealPtr stress, const FieldOnCellsRealPtr internVar,
    const ASTERDOUBLE &time_prev, const ASTERDOUBLE &time_step,
    const FieldOnCellsRealPtr &externVarPrev, const FieldOnCellsRealPtr &externVarCurr,
    const VectorString &groupOfCells ) const {
    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currBehaviour = _phys_problem->getBehaviourProperty();

    // Select option for matrix
    std::string option = "RIGI_MECA_TANG";

    // Prepare computing:
    CalculPtr calcul = createCalculForNonLinear( option, time_prev, time_prev + time_step,
                                                 externVarPrev, externVarCurr, groupOfCells );
    FiniteElementDescriptorPtr FEDesc = calcul->getFiniteElementDescriptor();

    // Set current physical state
    calcul->addInputField( "PDEPLMR", displ );
    calcul->addInputField( "PDEPLPR", displ_step );
    calcul->addInputField( "PCONTMR", stress );
    calcul->addInputField( "PVARIMR", internVar );

    // Provisoire: pour TANGENTE=VERIFICATION, nécessité de variables internes à chaque itération
    FieldOnCellsRealPtr vari_iter = std::make_shared< FieldOnCellsReal >(
        FEDesc, "ELGA", "VARI_R", currBehaviour, currElemChara );
    calcul->addInputField( "PVARIMP", vari_iter );

    // Create output matrix
    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    // Create output vector
    auto elemVect = std::make_shared< ElementaryVectorReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );
    elemVect->prepareCompute( option );

    // Create output fields
    FieldOnCellsRealPtr stress_pred =
        std::make_shared< FieldOnCellsReal >( currModel, "ELGA", "SIEF_R", currElemChara );
    FieldOnCellsLongPtr maskField = std::make_shared< FieldOnCellsLong >( FEDesc );
    FieldOnCellsLongPtr exitField = std::make_shared< FieldOnCellsLong >( FEDesc );

    // Add output fields
    calcul->addOutputField( "PCONTPR", stress_pred );
    calcul->addOutputField( "PCOPRED", maskField );
    calcul->addOutputField( "PCODRET", exitField );

    // Add output elementary
    calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
    calcul->addOutputElementaryTerm( "PMATUNS", std::make_shared< ElementaryTermReal >() );
    calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );

    // Compute
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUUR" ) );
        if ( calcul->hasOutputElementaryTerm( "PMATUNS" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUNS" ) );
        if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) )
            elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ) );
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

ElementaryMatrixDisplacementRealPtr DiscreteComputation::getContactMatrix(
    const MeshCoordinatesFieldPtr geom, const FieldOnNodesRealPtr displ,
    const FieldOnNodesRealPtr displ_step, const ASTERDOUBLE &time_prev,
    const ASTERDOUBLE &time_step, const FieldOnCellsRealPtr data,
    const FieldOnNodesRealPtr coef_cont, const FieldOnNodesRealPtr coef_frot ) const {
    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    // Select option for matrix
    std::string option = "RIGI_CONT";

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
    calcul->addInputField( "PCCONTR", coef_cont );
    calcul->addInputField( "PCFROTR", coef_frot );

    // Add time fields
    calcul->addTimeField( "PINSTMR", time_prev );
    calcul->addTimeField( "PINSTPR", time_prev + time_step );

    // Create output vector
    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    // Add output elementary
    calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
    calcul->addOutputElementaryTerm( "PMATUNS", std::make_shared< ElementaryTermReal >() );

    // Computation
    calcul->compute();
    if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) ) {
        elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUUR" ) );
    }
    if ( calcul->hasOutputElementaryTerm( "PMATUNS" ) ) {
        elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUNS" ) );
    }
    elemMatr->build();

    return elemMatr;
}

ElementaryMatrixDisplacementRealPtr
DiscreteComputation::getRotationalStiffnessMatrix( const VectorString &groupOfCells ) const {
    AS_ASSERT( _phys_problem->getModel()->isMechanical() );
    const std::string option = "RIGI_MECA_RO";

    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    // Prepare loads
    const auto listOfLoads = _phys_problem->getListOfLoads();
    const auto model = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Get gyroscopic fields
    std::vector< DataFieldPtr > rotat;
    auto impl = [&rotat]( auto loads ) {
        for ( const auto &load : loads ) {
            if ( load->hasLoadField( "ROTAT" ) ) {
                rotat.push_back( load->getConstantLoadField( "ROTAT" ) );
            }
        }
    };

    impl( listOfLoads->getMechanicalLoadsReal() );

    if ( rotat.empty() ) {
        rotat.push_back( nullptr );
    }

    if ( rotat.size() != 1 ) {
        UTMESS( "F", "CALCULEL3_71" );
    }

    for ( const auto &gyro : rotat ) {
        // Prepare computing
        auto calcul = std::make_unique< Calcul >( option );
        if ( groupOfCells.empty() ) {
            calcul->setModel( model );
        } else {
            calcul->setGroupsOfCells( model, groupOfCells );
        }

        // Add input fields
        calcul->addInputField( "PGEOMER", model->getMesh()->getCoordinates() );
        if ( currMater ) {
            calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
            calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );
        }
        if ( currElemChara ) {
            calcul->addElementaryCharacteristicsField( currElemChara );
        }

        if ( gyro ) {
            calcul->addInputField( "PROTATR", gyro );
        }

        calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) ) {
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUUR" ) );
        }
    }

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixDisplacementRealPtr
DiscreteComputation::getGyroscopicStiffnessMatrix( const VectorString &groupOfCells ) const {
    AS_ASSERT( _phys_problem->getModel()->isMechanical() );
    const std::string option = "RIGI_GYRO";

    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    // Prepare loads
    const auto listOfLoads = _phys_problem->getListOfLoads();
    const auto model = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Get gyroscopic fields
    std::vector< DataFieldPtr > rotat;
    auto impl = [&rotat]( auto loads ) {
        for ( const auto &load : loads ) {
            if ( load->hasLoadField( "ROTAT" ) ) {
                rotat.push_back( load->getConstantLoadField( "ROTAT" ) );
            }
        }
    };

    impl( listOfLoads->getMechanicalLoadsReal() );

    if ( rotat.empty() ) {
        rotat.push_back( nullptr );
    }

    if ( rotat.size() != 1 ) {
        UTMESS( "F", "CALCULEL3_71" );
    }

    for ( const auto &gyro : rotat ) {
        // Prepare computing
        auto calcul = std::make_unique< Calcul >( option );
        if ( groupOfCells.empty() ) {
            calcul->setModel( model );
        } else {
            calcul->setGroupsOfCells( model, groupOfCells );
        }

        // Add input fields
        calcul->addInputField( "PGEOMER", model->getMesh()->getCoordinates() );
        if ( currMater ) {
            calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
            calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );
        }
        if ( currElemChara ) {
            calcul->addElementaryCharacteristicsField( currElemChara );
        }

        if ( gyro ) {
            calcul->addInputField( "PROTATR", gyro );
        }

        calcul->addOutputElementaryTerm( "PMATUNS", std::make_shared< ElementaryTermReal >() );
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATUNS" ) ) {
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUNS" ) );
        }
    }

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixDisplacementRealPtr
DiscreteComputation::getGyroscopicDampingMatrix( const VectorString &groupOfCells ) const {
    AS_ASSERT( _phys_problem->getModel()->isMechanical() );
    const std::string option = "MECA_GYRO";

    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    // Prepare loads
    const auto listOfLoads = _phys_problem->getListOfLoads();
    const auto model = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Get gyroscopic fields
    std::vector< DataFieldPtr > rotat;
    auto impl = [&rotat]( auto loads ) {
        for ( const auto &load : loads ) {
            if ( load->hasLoadField( "ROTAT" ) ) {
                rotat.push_back( load->getConstantLoadField( "ROTAT" ) );
            }
        }
    };

    impl( listOfLoads->getMechanicalLoadsReal() );

    if ( rotat.empty() ) {
        rotat.push_back( nullptr );
    }

    if ( rotat.size() != 1 ) {
        UTMESS( "F", "CALCULEL3_71" );
    }

    for ( const auto &gyro : rotat ) {
        // Prepare computing
        auto calcul = std::make_unique< Calcul >( option );
        if ( groupOfCells.empty() ) {
            calcul->setModel( model );
        } else {
            calcul->setGroupsOfCells( model, groupOfCells );
        }

        // Add input fields
        calcul->addInputField( "PGEOMER", model->getMesh()->getCoordinates() );
        if ( currMater ) {
            calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
            calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );
        }

        if ( currElemChara ) {
            calcul->addElementaryCharacteristicsField( currElemChara );
        }

        if ( gyro ) {
            calcul->addInputField( "PROTATR", gyro );
        }

        calcul->addOutputElementaryTerm( "PMATUNS", std::make_shared< ElementaryTermReal >() );
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATUNS" ) ) {
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUNS" ) );
        }
    }

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixDisplacementRealPtr
DiscreteComputation::getImpedanceBoundaryMatrix( const VectorString &groupOfCells,
                                                 const ASTERINTEGER &onde_flui ) const {
    AS_ASSERT( _phys_problem->getModel()->isMechanical() );
    const std::string option = "IMPE_MECA";

    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    const auto model = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();

    auto calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( model );
    } else {
        calcul->setGroupsOfCells( model, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", model->getMesh()->getCoordinates() );
    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
    }
    auto carteFluid = createWaveTypeFluidField( onde_flui );
    calcul->addInputField( "PWAVETFL", carteFluid );

    calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
    calcul->compute();
    if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) ) {
        elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUUR" ) );
    }

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixDisplacementRealPtr
DiscreteComputation::getImpedanceWaveMatrix( const VectorString &groupOfCells ) const {
    AS_ASSERT( _phys_problem->getModel()->isMechanical() );
    const std::string option = "ONDE_FLUI";

    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics() );
    elemMatr->prepareCompute( option );

    // Prepare loads
    const auto listOfLoads = _phys_problem->getListOfLoads();
    const auto model = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();

    // Get wave fields
    std::vector< std::pair< std::string, DataFieldPtr > > wave;
    auto impl = [&wave]( auto loads, std::string param ) {
        for ( const auto &load : loads ) {
            if ( load->hasLoadField( "ONDE" ) ) {
                wave.push_back( std::make_pair( param, load->getConstantLoadField( "ONDE" ) ) );
            }
        }
    };

    impl( listOfLoads->getMechanicalLoadsReal(), "PONDECR" );

    if ( wave.empty() ) {
        UTMESS( "F", "CHARGES6_84" );
    }

    for ( const auto &[param, field] : wave ) {
        // Prepare computing
        auto calcul = std::make_unique< Calcul >( option );
        if ( groupOfCells.empty() ) {
            calcul->setModel( model );
        } else {
            calcul->setGroupsOfCells( model, groupOfCells );
        }

        // Add input fields
        calcul->addInputField( "PGEOMER", model->getMesh()->getCoordinates() );
        if ( currMater ) {
            calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
        }

        calcul->addInputField( param, field );

        calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATUUR" ) ) {
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATUUR" ) );
        }
    }

    elemMatr->build();
    return elemMatr;
};
