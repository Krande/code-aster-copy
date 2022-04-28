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

FieldOnNodesRealPtr DiscreteComputation::imposedDisplacement( ASTERDOUBLE currTime ) {

    ElementaryVectorDisplacementRealPtr elemVect =
        std::make_shared< ElementaryVectorDisplacementReal >();

    if ( _phys_problem->getModel()->isThermal() ) {
        AS_ABORT( "Not implemented for thermic" );
    };

    // Prepare loads
    auto listOfLoads = _phys_problem->getListOfLoads();
    JeveuxVectorChar24 listOfLoadsList = listOfLoads->getListVector();
    JeveuxVectorLong listOfLoadsInfo = listOfLoads->getInformationVector();

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( _phys_problem->getModel()->getName(), 24 );
    std::string nameLcha = ljust( listOfLoadsList->getName(), 24 );
    std::string nameInfc = ljust( listOfLoadsInfo->getName(), 24 );
    std::string vectElemName = ljust( elemVect->getName(), 24 );
    std::string typres( "R" );

    // Wrapper FORTRAN
    CALLO_VEDIME( modelName, nameLcha, nameInfc, &currTime, typres, vectElemName );

    // Construct vect_elem object
    elemVect->setListOfLoads( listOfLoads );
    elemVect->setModel( _phys_problem->getModel() );
    elemVect->build();

    // Assemble
    return elemVect->assembleWithLoadFunctions( _phys_problem->getDOFNumbering(), currTime );
};

FieldOnNodesRealPtr DiscreteComputation::dualReaction( FieldOnNodesRealPtr lagr_curr ) {

    ElementaryVectorDisplacementRealPtr elemVect =
        std::make_shared< ElementaryVectorDisplacementReal >();

    if ( _phys_problem->getModel()->isThermal() ) {
        AS_ABORT( "Not implemented for thermic" );
    };

    // Prepare loads
    auto listOfLoads = _phys_problem->getListOfLoads();
    std::string listLoadsName = ljust( listOfLoads->getName(), 19 );

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( _phys_problem->getModel()->getName(), 24 );
    std::string materName = ljust( _phys_problem->getMaterialField()->getName(), 24 );
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    std::string caraName( " " );
    if ( currElemChara )
        caraName = currElemChara->getName();
    caraName.resize( 24, ' ' );
    std::string vectElemName = ljust( elemVect->getName(), 24 );
    std::string base( "G" );
    std::string lagrName = lagr_curr->getName();

    // Wrapper FORTRAN
    CALLO_VEBTLA( base, modelName, materName, caraName, lagrName, listLoadsName, vectElemName );

    // Construct vect_elem object
    elemVect->setListOfLoads( listOfLoads );
    elemVect->setModel( _phys_problem->getModel() );
    elemVect->build();

    // Assemble
    return elemVect->assemble( _phys_problem->getDOFNumbering() );
};

FieldOnNodesRealPtr DiscreteComputation::dualDisplacement( FieldOnNodesRealPtr disp_curr,
                                                           ASTERDOUBLE scaling ) {

    ElementaryVectorDisplacementRealPtr elemVect =
        std::make_shared< ElementaryVectorDisplacementReal >();

    if ( _phys_problem->getModel()->isThermal() ) {
        AS_ABORT( "Not implemented for thermic" );
    };

    // Prepare loads
    auto listOfLoads = _phys_problem->getListOfLoads();
    std::string listLoadsName = ljust( listOfLoads->getName(), 19 );

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( _phys_problem->getModel()->getName(), 24 );
    std::string dispName = ljust( disp_curr->getName(), 24 );
    std::string vectElemName = ljust( elemVect->getName(), 24 );
    const std::string base( "G" );
    const ASTERDOUBLE const_scaling = scaling;

    // Wrapper FORTRAN
    CALLO_VEBUME( modelName, dispName, listLoadsName, vectElemName, &const_scaling, base );

    // Construct vect_elem object
    elemVect->setListOfLoads( listOfLoads );
    elemVect->setModel( _phys_problem->getModel() );
    elemVect->build();

    // Assemble
    FieldOnNodesRealPtr bume = elemVect->assemble( _phys_problem->getDOFNumbering() );

    if ( _phys_problem->getMesh()->isParallel() )
        CALLO_AP_ASSEMBLY_VECTOR( bume->getName() );

    return bume;
};

FieldOnNodesRealPtr DiscreteComputation::neumann( const VectorReal timeParameters,
                                                  const FieldOnCellsRealPtr _externVarField ) {

    ElementaryVectorDisplacementRealPtr elemVect =
        std::make_shared< ElementaryVectorDisplacementReal >();

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto listOfLoads = _phys_problem->getListOfLoads();

    if ( currModel->isThermal() ) {
        AS_ABORT( "Not implemented for thermic" );
    };

    if ( timeParameters.size() != 3 )
        throw std::runtime_error( "Invalid number of parameter" );
    const ASTERDOUBLE &currTime = timeParameters[0];

    // Prepare loads
    JeveuxVectorChar24 listOfLoadsList = listOfLoads->getListVector();
    JeveuxVectorLong listOfLoadsInfo = listOfLoads->getInformationVector();
    std::string nameLcha = ljust( listOfLoadsList->getName(), 24 );
    std::string nameInfc = ljust( listOfLoadsInfo->getName(), 24 );

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( currModel->getName(), 24 );
    std::string materName = ljust( currMater->getName(), 24 );
    std::string currCodedMaterName = ljust( currCodedMater->getName(), 24 );
    std::string currElemCharaName( " " );
    if ( currElemChara )
        currElemCharaName = currElemChara->getName();
    currElemCharaName.resize( 24, ' ' );
    std::string vectElemName = ljust( elemVect->getName(), 24 );
    std::string stop( "S" );

    // Get external state variables
    std::string externVarName( " " );
    if ( _externVarField ) {
        externVarName = _externVarField->getName();
    }
    externVarName.resize( 24, ' ' );

    // Wrapper FORTRAN
    CALLO_VECHME_WRAP( stop, modelName, nameLcha, nameInfc, &currTime, currElemCharaName, materName,
                       currCodedMaterName, vectElemName, externVarName );

    // Construct vect_elem object
    elemVect->setListOfLoads( listOfLoads );
    elemVect->setModel( _phys_problem->getModel() );
    elemVect->build();

    // Assemble
    return elemVect->assembleWithLoadFunctions( _phys_problem->getDOFNumbering(),
                                                timeParameters[0] + timeParameters[1] );
};

FieldOnNodesRealPtr DiscreteComputation::dirichletBC( const ASTERDOUBLE &time ) const {

    auto dofNume = _phys_problem->getDOFNumbering();
    FieldOnNodesRealPtr vectAsse = std::make_shared< FieldOnNodesReal >( dofNume );

    // Prepare loads
    const auto &_listOfLoads = _phys_problem->getListOfLoads();
    if ( _listOfLoads->isEmpty() )
        _listOfLoads->build( _phys_problem->getModel() );

    JeveuxVectorChar24 listOfLoadsList = _listOfLoads->getListVector();
    JeveuxVectorLong listOfLoadsInfo = _listOfLoads->getInformationVector();
    JeveuxVectorChar24 listOfLoadsFunc = _listOfLoads->getListOfFunctions();
    std::string nameLcha = ljust( listOfLoadsList->getName(), 24 );
    std::string nameInfc = ljust( listOfLoadsInfo->getName(), 24 );
    std::string nameFcha = ljust( listOfLoadsFunc->getName(), 24 );

    // Get JEVEUX names of objects to call Fortran
    std::string vectAsseName = vectAsse->getName();
    std::string dofNumName = _phys_problem->getDOFNumbering()->getName();
    std::string base( "G" );

    // Wrapper FORTRAN
    CALLO_ASCAVC_WRAP( nameLcha, nameInfc, nameFcha, dofNumName, &time, vectAsseName, base );

    // Construct vect_asse object
    vectAsse->build();

    // Assemble
    return vectAsse;
};

FieldOnNodesRealPtr
DiscreteComputation::incrementalDirichletBC( const ASTERDOUBLE &time,
                                             const FieldOnNodesRealPtr disp_curr ) const {
    auto dofNume = _phys_problem->getDOFNumbering();

    if ( dofNume->hasDirichletBC() ) {
        auto diri_curr = dirichletBC( time );
        auto diri_impo = *( diri_curr ) - *( disp_curr );

        // Set to zero terms not imposed
        auto eliminatedDofs = dofNume->getDirichletBCDOFs();
        auto nbElimination = eliminatedDofs.size();

        for ( ASTERINTEGER ieq = 0; ieq < nbElimination; ieq++ ) {
            if ( eliminatedDofs[ieq] == 0 )
                diri_impo[ieq] = 0.0;
        }

        return std::make_shared< FieldOnNodesReal >( diri_impo );
    }

    // Construct vect_asse object
    FieldOnNodesRealPtr vectAsse = std::make_shared< FieldOnNodesReal >( dofNume );
    vectAsse->setValues( 0.0 );
    vectAsse->build();

    return vectAsse;
};

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

ElementaryMatrixDisplacementRealPtr DiscreteComputation::massMatrix( ASTERDOUBLE time ) {
    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >();

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

FieldOnNodesRealPtr DiscreteComputation::computeExternalStateVariablesLoad(
    const ASTERDOUBLE &time, const ConstantFieldOnCellsRealPtr _timeField,
    const FieldOnCellsRealPtr _externVarField ) const {

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currExternVarRefe = _phys_problem->getExternalStateVariablesReference();

    // Some checks
    AS_ASSERT( currMater );
    AS_ASSERT( currMater->hasExternalStateVariableForLoad() );
    if ( currMater->hasExternalStateVariableWithReference() ) {
        AS_ASSERT( currExternVarRefe );
    }

    // Main object
    CalculPtr _calcul = std::make_unique< Calcul >( "CHAR_VARC" );

    // Create specific output field for XFEM
    FieldOnCellsRealPtr sigmXfem;
    if ( currModel->existsXfem() ) {
        const std::string option = "SIEF_ELGA";
        const std::string paraName = "PCONTRR";
        sigmXfem = std::make_shared< FieldOnCellsReal >( currModel, option, paraName );
    }

    // Create elementary vectors
    ElementaryVectorDisplacementRealPtr elemVect =
        std::make_shared< ElementaryVectorDisplacementReal >();
    elemVect->setModel( currModel );
    elemVect->setMaterialField( currMater );
    if ( currElemChara ) {
        elemVect->setElementaryCharacteristics( currElemChara );
    }
    elemVect->prepareCompute( "CHAR_VARC" );
    int nbExternVar = static_cast< int >( externVarEnumInt::NumberOfExternVarTypes );
    for ( auto iExternVar = 0; iExternVar < nbExternVar; iExternVar++ ) {
        externVarEnumInt numeExternVar = static_cast< externVarEnumInt >( iExternVar );
        if ( currMater->hasExternalStateVariable( numeExternVar ) &&
             ExternalVariableTraits::externVarHasStrain( numeExternVar ) ) {
            const auto option = ExternalVariableTraits::getExternVarOption( numeExternVar );
            _calcul->setOption( option );
            _calcul->setModel( currModel );
            _calcul->clearInputs();
            _calcul->clearOutputs();

            // Add input fields
            _calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            _calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
            _calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );
            if ( currElemChara ) {
                _calcul->addElementaryCharacteristicsField( currElemChara );
            }
            _calcul->addInputField( "PVARCPR", _externVarField );
            _calcul->addInputField( "PTEMPSR", _timeField );
            if ( currExternVarRefe ) {
                _calcul->addInputField( "PVARCRR", currExternVarRefe );
            }
            if ( currModel->existsXfem() ) {
                XfemModelPtr currXfemModel = currModel->getXfemModel();
                _calcul->addXFEMField( currXfemModel );
                _calcul->addOutputField( "PCONTRT", sigmXfem );
            }
            _calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );
            _calcul->compute();
            if ( _calcul->hasOutputElementaryTerm( "PVECTUR" ) ) {
                elemVect->addElementaryTerm( _calcul->getOutputElementaryTerm( "PVECTUR" ) );
            }
        }
    }
    // Build elementary vectors
    elemVect->build();

    // Assemble
    return elemVect->assemble( _phys_problem->getDOFNumbering() );
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
    const FieldOnCellsRealPtr _externVarFieldCurr ) {

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
    _calcul->setModel( currModel );

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

/** @brief Compute internal forces, stress and internal state variables */
std::tuple< FieldOnCellsLongPtr, FieldOnCellsRealPtr, FieldOnCellsRealPtr,
            ElementaryVectorDisplacementRealPtr >
DiscreteComputation::computeInternalForces( const FieldOnNodesRealPtr displ,
                                            const FieldOnNodesRealPtr displ_incr,
                                            const FieldOnCellsRealPtr stress,
                                            const FieldOnCellsRealPtr _internVar,
                                            const ConstantFieldOnCellsRealPtr _timeFieldPrev,
                                            const ConstantFieldOnCellsRealPtr _timeFieldCurr ) {

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
    CalculPtr _calcul = createCalculForNonLinear( option, _timeFieldPrev, _timeFieldCurr,
                                                  _externVarFieldPrev, _externVarFieldCurr );

    // Set current physical state
    _calcul->addInputField( "PDEPLMR", displ );
    _calcul->addInputField( "PDEPLPR", displ_incr );
    _calcul->addInputField( "PCONTMR", stress );
    _calcul->addInputField( "PVARIMR", _internVar );

    // Provisoire: pour TANGENTE=VERIFICATION, nécessité de variables internes à chaque itération
    FieldOnCellsRealPtr vari_iter = std::make_shared< FieldOnCellsReal >(
        currModel, currBehaviour, "ELGA_VARI_R", currElemChara );
    _calcul->addInputField( "PVARIMP", vari_iter );

    // Create output vector
    auto elemVect = std::make_shared< ElementaryVectorDisplacementReal >();
    elemVect->setModel( currModel );
    elemVect->setMaterialField( currMater );
    elemVect->setElementaryCharacteristics( currElemChara );
    elemVect->prepareCompute( option );

    // Create output fields
    FieldOnCellsRealPtr stress_curr =
        std::make_shared< FieldOnCellsReal >( currModel, nullptr, "ELGA_SIEF_R", currElemChara );
    FieldOnCellsLongPtr exitField = std::make_shared< FieldOnCellsLong >();
    FieldOnCellsRealPtr vari_curr = std::make_shared< FieldOnCellsReal >(
        currModel, currBehaviour, "ELGA_VARI_R", currElemChara );

    // Add output fields
    _calcul->addOutputField( "PVARIPR", vari_curr );
    _calcul->addOutputField( "PCONTPR", stress_curr );
    _calcul->addOutputField( "PCODRET", exitField );

    // Add output elementary
    _calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );

    // Compute
    if ( currModel->existsFiniteElement() ) {
        _calcul->compute();
        if ( _calcul->hasOutputElementaryTerm( "PVECTUR" ) )
            elemVect->addElementaryTerm( _calcul->getOutputElementaryTerm( "PVECTUR" ) );
        elemVect->build();
    };

    return std::make_tuple( exitField, vari_curr, stress_curr, elemVect );
}

/** @brief Compute tangent matrix (not assembled) */
std::tuple< FieldOnCellsLongPtr, FieldOnCellsRealPtr, FieldOnCellsRealPtr,
            ElementaryVectorDisplacementRealPtr, ElementaryMatrixDisplacementRealPtr >
DiscreteComputation::computeTangentStiffnessMatrix(
    const FieldOnNodesRealPtr displ, const FieldOnNodesRealPtr displ_incr,
    const FieldOnCellsRealPtr stress, const FieldOnCellsRealPtr _internVar,
    const ConstantFieldOnCellsRealPtr _timeFieldPrev,
    const ConstantFieldOnCellsRealPtr _timeFieldCurr ) {

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
    CalculPtr _calcul = createCalculForNonLinear( option, _timeFieldPrev, _timeFieldCurr,
                                                  _externVarFieldPrev, _externVarFieldCurr );

    // Set current physical state
    _calcul->addInputField( "PDEPLMR", displ );
    _calcul->addInputField( "PDEPLPR", displ_incr );
    _calcul->addInputField( "PCONTMR", stress );
    _calcul->addInputField( "PVARIMR", _internVar );

    // Provisoire: pour TANGENTE=VERIFICATION, nécessité de variables internes à chaque itération
    FieldOnCellsRealPtr vari_iter = std::make_shared< FieldOnCellsReal >(
        currModel, currBehaviour, "ELGA_VARI_R", currElemChara );
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
    FieldOnCellsRealPtr stress_curr =
        std::make_shared< FieldOnCellsReal >( currModel, nullptr, "ELGA_SIEF_R", currElemChara );
    FieldOnCellsLongPtr exitField = std::make_shared< FieldOnCellsLong >();
    FieldOnCellsRealPtr vari_curr = std::make_shared< FieldOnCellsReal >(
        currModel, currBehaviour, "ELGA_VARI_R", currElemChara );

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
    return std::make_tuple( exitField, vari_curr, stress_curr, elemVect, elemMatr );
}

/** @brief Compute tangent prediction matrix (not assembled) */
std::tuple< FieldOnCellsLongPtr, FieldOnCellsRealPtr, ElementaryMatrixDisplacementRealPtr,
            ElementaryVectorDisplacementRealPtr >
DiscreteComputation::computeTangentPredictionMatrix(
    const FieldOnNodesRealPtr displ, const FieldOnNodesRealPtr displ_incr,
    const FieldOnCellsRealPtr stress, const FieldOnCellsRealPtr _internVar,
    const ConstantFieldOnCellsRealPtr _timeFieldPrev,
    const ConstantFieldOnCellsRealPtr _timeFieldCurr ) {

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
    CalculPtr _calcul = createCalculForNonLinear( option, _timeFieldPrev, _timeFieldCurr,
                                                  _externVarFieldPrev, _externVarFieldCurr );

    // Set current physical state
    _calcul->addInputField( "PDEPLMR", displ );
    _calcul->addInputField( "PDEPLPR", displ_incr );
    _calcul->addInputField( "PCONTMR", stress );
    _calcul->addInputField( "PVARIMR", _internVar );

    // Provisoire: pour TANGENTE=VERIFICATION, nécessité de variables internes à chaque itération
    FieldOnCellsRealPtr vari_iter = std::make_shared< FieldOnCellsReal >(
        currModel, currBehaviour, "ELGA_VARI_R", currElemChara );
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
    FieldOnCellsRealPtr stress_pred =
        std::make_shared< FieldOnCellsReal >( currModel, nullptr, "ELGA_SIEF_R", currElemChara );
    FieldOnCellsLongPtr maskField = std::make_shared< FieldOnCellsLong >();
    FieldOnCellsLongPtr exitField = std::make_shared< FieldOnCellsLong >();

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
    return std::make_tuple( maskField, stress_pred, elemMatr, elemVect );
}
