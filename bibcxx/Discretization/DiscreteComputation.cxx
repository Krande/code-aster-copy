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

#include "Discretization/DiscreteComputation.h"

#include "aster_fort_calcul.h"
#include "aster_fort_superv.h"
#include "astercxx.h"

#include "Discretization/Calcul.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Modeling/XfemModel.h"
#include "Utilities/Tools.h"


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
ConstantFieldOnCellsRealPtr DiscreteComputation::createTimeField( const ASTERDOUBLE time ) {

    // Get mesh
    auto mesh = _phys_problem->getMesh();

    // Create field
    auto field = std::make_shared< ConstantFieldOnCellsReal >( mesh );

    // Get JEVEUX names of objects to call Fortran
    const std::string physicalName( "INST_R" );
    field->allocate( physicalName );
    ConstantFieldOnZone a( mesh );
    ConstantFieldValues< ASTERDOUBLE > b( { "INST" }, { time } );
    field->setValueOnZone( a, b );

    return field;
}

FieldOnCellsRealPtr
DiscreteComputation::createExternalStateVariablesField( const ASTERDOUBLE time ) {

    // Create field
    auto field = std::make_shared< FieldOnCellsReal >();

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( _phys_problem->getModel()->getName(), 24 );
    std::string materialFieldName = ljust( _phys_problem->getMaterialField()->getName(), 24 );
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    std::string elemCharaName( " " );
    if ( currElemChara )
        elemCharaName = currElemChara->getName();
    elemCharaName.resize( 24, ' ' );
    std::string fieldName = ljust( field->getName(), 19 );

    // Output
    std::string out( ' ', 2 );

    // Call Fortran WRAPPER
    CALLO_VRCINS_WRAP( modelName, materialFieldName, elemCharaName, &time, fieldName, out );

    return field;
}

FieldOnCellsRealPtr DiscreteComputation::computeExternalStateVariablesReference() const {

    // Create field
    auto field = std::make_shared< FieldOnCellsReal >();

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( _phys_problem->getModel()->getName(), 8 );
    std::string materialFieldName = ljust( _phys_problem->getMaterialField()->getName(), 8 );
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    std::string elemCharaName( " ", 8 );
    if ( currElemChara )
        elemCharaName = std::string( currElemChara->getName(), 0, 8 );
    std::string fieldName = ljust( field->getName(), 19 );

    // Call Fortran WRAPPER
    CALLO_VRCREF( modelName, materialFieldName, elemCharaName, fieldName );

    return field;
}

CalculPtr DiscreteComputation::createCalculForNonLinear(
    const std::string option, const ConstantFieldOnCellsRealPtr _timeFieldPrev,
    const ConstantFieldOnCellsRealPtr _timeFieldCurr, const FieldOnCellsRealPtr _externVarFieldPrev,
    const FieldOnCellsRealPtr _externVarFieldCurr, const VectorString &groupOfCells ) {

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currBehaviour = _phys_problem->getBehaviourProperty();
    auto currExternVarRefe = _phys_problem->getExternalStateVariablesReference();
    AS_ASSERT( currMater );

    // No !
    if ( currModel->existsHHO() ) {
        throw std::runtime_error( "HHO not implemented" );
    }
    if ( currModel->exists3DShell() ) {
        throw std::runtime_error( "COQUE_3D not implemented" );
    }
    if ( currModel->existsSTRX() ) {
        throw std::runtime_error( "Beams not implemented" );
    }

    // Prepare computing: the main object
    CalculPtr _calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        _calcul->setModel( currModel );
    } else {
        _calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add external state variables
    if ( currMater->hasExternalStateVariableWithReference() ) {
        AS_ASSERT( currExternVarRefe );
        _calcul->addInputField( "PVARCRR", currExternVarRefe );
    }
    if ( currMater->hasExternalStateVariable() ) {
        if ( !_externVarFieldPrev ) {
            AS_ABORT( "External state variables vector for beginning of time step is missing" )
        }
        if ( !_externVarFieldCurr ) {
            AS_ABORT( "External state variables vector for end of time step is missing" )
        }
        AS_ASSERT( currExternVarRefe );
        _calcul->addInputField( "PVARCMR", _externVarFieldPrev );
        _calcul->addInputField( "PVARCPR", _externVarFieldCurr );
    }

    // Add time fields
    if ( !_timeFieldPrev ) {
        AS_ABORT( "Time field for beginning of time step is missing" )
    }
    if ( !_timeFieldCurr ) {
        AS_ABORT( "Time field for end of time step is missing" )
    }
    _calcul->addInputField( "PINSTMR", _timeFieldPrev );
    _calcul->addInputField( "PINSTPR", _timeFieldCurr );

    // Add input fields
    _calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    _calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
    if ( currElemChara ) {
        _calcul->addElementaryCharacteristicsField( currElemChara );
    }
    if ( currModel->existsXfem() ) {
        XfemModelPtr currXfemModel = currModel->getXfemModel();
        _calcul->addXFEMField( currXfemModel );
    }
    _calcul->addBehaviourField( currBehaviour );

    return _calcul;
};
