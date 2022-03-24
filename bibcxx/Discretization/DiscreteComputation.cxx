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

    if ( _study->getModel()->isThermal() ) {
        AS_ABORT( "Not implemented for thermic" );
    };

    // Prepare loads
    auto listOfLoads = _study->getListOfLoads();
    JeveuxVectorChar24 listOfLoadsList = listOfLoads->getListVector();
    JeveuxVectorLong listOfLoadsInfo = listOfLoads->getInformationVector();

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( _study->getModel()->getName(), 24 );
    std::string nameLcha = ljust( listOfLoadsList->getName(), 24 );
    std::string nameInfc = ljust( listOfLoadsInfo->getName(), 24 );
    std::string vectElemName = ljust( elemVect->getName(), 24 );
    std::string typres( "R" );

    // Wrapper FORTRAN
    CALLO_VEDIME( modelName, nameLcha, nameInfc, &currTime, typres, vectElemName );

    // Construct vect_elem object
    elemVect->setListOfLoads( listOfLoads );
    elemVect->setModel( _study->getModel() );
    elemVect->build();

    // Assemble
    return elemVect->assembleWithLoadFunctions( _study->getDOFNumbering(), currTime );
};

FieldOnNodesRealPtr DiscreteComputation::dualReaction( FieldOnNodesRealPtr lagr_curr ) {

    ElementaryVectorDisplacementRealPtr elemVect =
        std::make_shared< ElementaryVectorDisplacementReal >();

    if ( _study->getModel()->isThermal() ) {
        AS_ABORT( "Not implemented for thermic" );
    };

    // Prepare loads
    auto listOfLoads = _study->getListOfLoads();
    std::string listLoadsName = ljust( listOfLoads->getName(), 19 );

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( _study->getModel()->getName(), 24 );
    std::string materName = ljust( _study->getMaterialField()->getName(), 24 );
    auto currElemChara = _study->getElementaryCharacteristics();
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
    elemVect->setModel( _study->getModel() );
    elemVect->build();

    // Assemble
    return elemVect->assemble( _study->getDOFNumbering() );
};

FieldOnNodesRealPtr DiscreteComputation::dualDisplacement( FieldOnNodesRealPtr disp_curr,
                                                           ASTERDOUBLE scaling ) {

    ElementaryVectorDisplacementRealPtr elemVect =
        std::make_shared< ElementaryVectorDisplacementReal >();

    if ( _study->getModel()->isThermal() ) {
        AS_ABORT( "Not implemented for thermic" );
    };

    // Prepare loads
    auto listOfLoads = _study->getListOfLoads();
    std::string listLoadsName = ljust( listOfLoads->getName(), 19 );

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( _study->getModel()->getName(), 24 );
    std::string dispName = ljust( disp_curr->getName(), 24 );
    std::string vectElemName = ljust( elemVect->getName(), 24 );
    const std::string base( "G" );
    const ASTERDOUBLE const_scaling = scaling;

    // Wrapper FORTRAN
    CALLO_VEBUME( modelName, dispName, listLoadsName, vectElemName, &const_scaling, base );

    // Construct vect_elem object
    elemVect->setListOfLoads( listOfLoads );
    elemVect->setModel( _study->getModel() );
    elemVect->build();

    // Assemble
    FieldOnNodesRealPtr bume = elemVect->assemble( _study->getDOFNumbering() );

    if ( _study->getMesh()->isParallel() )
        CALLO_AP_ASSEMBLY_VECTOR( bume->getName() );

    return bume;
};

FieldOnNodesRealPtr DiscreteComputation::neumann( const VectorReal timeParameters,
                                                  const FieldOnCellsRealPtr _externVarField ) {

    ElementaryVectorDisplacementRealPtr elemVect =
        std::make_shared< ElementaryVectorDisplacementReal >();

    // Get main parameters
    auto currModel = _study->getModel();
    auto currMater = _study->getMaterialField();
    auto currCodedMater = _study->getCodedMaterial();
    auto currElemChara = _study->getElementaryCharacteristics();
    auto listOfLoads = _study->getListOfLoads();

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
    elemVect->setModel( _study->getModel() );
    elemVect->build();

    // Assemble
    return elemVect->assembleWithLoadFunctions( _study->getDOFNumbering(),
                                                timeParameters[0] + timeParameters[1] );
};

FieldOnNodesRealPtr DiscreteComputation::dirichletBC( const ASTERDOUBLE &time ) const {

    auto dofNume = _study->getDOFNumbering();
    FieldOnNodesRealPtr vectAsse = std::make_shared< FieldOnNodesReal >( dofNume );

    // Prepare loads
    const auto &_listOfLoads = _study->getListOfLoads();
    if ( _listOfLoads->isEmpty() )
        _listOfLoads->build( _study->getModel() );

    JeveuxVectorChar24 listOfLoadsList = _listOfLoads->getListVector();
    JeveuxVectorLong listOfLoadsInfo = _listOfLoads->getInformationVector();
    JeveuxVectorChar24 listOfLoadsFunc = _listOfLoads->getListOfFunctions();
    std::string nameLcha = ljust( listOfLoadsList->getName(), 24 );
    std::string nameInfc = ljust( listOfLoadsInfo->getName(), 24 );
    std::string nameFcha = ljust( listOfLoadsFunc->getName(), 24 );

    // Get JEVEUX names of objects to call Fortran
    std::string vectAsseName = vectAsse->getName();
    std::string dofNumName = _study->getDOFNumbering()->getName();
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
    auto dofNume = _study->getDOFNumbering();

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
    auto currModel = _study->getModel();
    auto currMater = _study->getMaterialField();
    auto currCodedMater = _study->getCodedMaterial();
    auto currElemChara = _study->getElementaryCharacteristics();

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
    _calcul->addTimeField( time );
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
    const auto &_listOfLoads = _study->getListOfLoads();

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
    ModelPtr currModel = _study->getModel();
    MaterialFieldPtr currMater = _study->getMaterialField();
    ElementaryCharacteristicsPtr currElemChara = _study->getElementaryCharacteristics();

    // Set parameters of elementary matrix
    elemMatr->setModel( currModel );
    elemMatr->setMaterialField( currMater );
    elemMatr->setElementaryCharacteristics( currElemChara );

    // Prepare computing
    const std::string option( "MECA_DDLM_R" );
    CalculPtr _calcul = std::make_shared< Calcul >( option );
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

ConstantFieldOnCellsRealPtr DiscreteComputation::createTimeField( const std::string fieldName,
                                                                  const ASTERDOUBLE time ) {

    // Get mesh
    auto mesh = _study->getMesh();

    // Create field
    auto field = std::make_shared< ConstantFieldOnCellsReal >( fieldName, mesh );

    // Get JEVEUX names of objects to call Fortran
    const std::string physicalName( "INST_R" );
    field->allocate( physicalName );
    ConstantFieldOnZone a( mesh );
    ConstantFieldValues< ASTERDOUBLE > b( { "INST" }, { time } );
    field->setValueOnZone( a, b );

    return field;
}

FieldOnCellsRealPtr
DiscreteComputation::createExternalStateVariablesField( const std::string fieldName,
                                                        const ASTERDOUBLE time ) {

    // Create field
    auto field = std::make_shared< FieldOnCellsReal >( fieldName );

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( _study->getModel()->getName(), 24 );
    std::string materialFieldName = ljust( _study->getMaterialField()->getName(), 24 );
    auto currElemChara = _study->getElementaryCharacteristics();
    std::string elemCharaName( " " );
    if ( currElemChara )
        elemCharaName = currElemChara->getName();
    elemCharaName.resize( 24, ' ' );

    // Output
    std::string out( ' ', 2 );

    // Call Fortran WRAPPER
    CALLO_VRCINS_WRAP( modelName, materialFieldName, elemCharaName, &time, fieldName, out );

    return field;
}

FieldOnNodesRealPtr DiscreteComputation::computeExternalStateVariablesLoad() const {

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( _study->getModel()->getName(), 24 );
    std::string materialFieldName = ljust( _study->getMaterialField()->getName(), 8 );
    std::string behaviourName =
        ljust( _study->getMaterialField()->getBehaviourField()->getName(), 24 );
    const auto compor = _study->getMaterialField()->getBehaviourField();
    std::string codedMaterialFieldName =
        ljust( _study->getCodedMaterial()->getCodedMaterialField()->getName(), 24 );
    auto currElemChara = _study->getElementaryCharacteristics();
    std::string elemCharaName( " " );
    if ( currElemChara )
        elemCharaName = currElemChara->getName();
    elemCharaName.resize( 24, ' ' );
    auto dofNumbering = _study->getDOFNumbering();
    std::string dofNumberingName( dofNumbering->getName() );

    // Detect external state variables
    ASTERINTEGER hydr = 0, sech = 0, temp = 0, ptot = 0;
    if ( _study->getMaterialField()->hasExternalStateVariable( externVarEnumInt::Temperature ) )
        temp = 1;
    if ( _study->getMaterialField()->hasExternalStateVariable(
             externVarEnumInt::ConcreteHydratation ) )
        hydr = 1;
    if ( _study->getMaterialField()->hasExternalStateVariable( externVarEnumInt::ConcreteDrying ) )
        sech = 1;
    if ( _study->getMaterialField()->hasExternalStateVariable(
             externVarEnumInt::TotalFluidPressure ) )
        ptot = 1;

    // Call Fortran WRAPPER (oui, c'est pourri, le nom est en dur)
    auto exteStateVariName = materialFieldName;
    std::string out( 24, ' ' );
    CALLO_CACHVC( modelName, materialFieldName, codedMaterialFieldName, elemCharaName,
                  dofNumberingName, behaviourName, exteStateVariName, out, &hydr, &sech, &temp,
                  &ptot );

    // Get name of nodal field from Fortran
    JeveuxVectorChar24 vectOut( out );
    vectOut->updateValuePointer();

    // Create nodal field
    FieldOnNodesRealPtr exteVariLoad( new FieldOnNodesReal( ( *vectOut )[0].toString() ) );
    exteVariLoad->setDescription( dofNumbering->getDescription() );
    exteVariLoad->setMesh( dofNumbering->getMesh() );
    exteVariLoad->build();
    exteVariLoad->updateValuePointers();

    return exteVariLoad;
}

FieldOnCellsRealPtr
DiscreteComputation::computeExternalStateVariablesReference( const std::string fieldName ) const {

    // Create field
    auto field = std::make_shared< FieldOnCellsReal >( fieldName );

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( _study->getModel()->getName(), 8 );
    std::string materialFieldName = ljust( _study->getMaterialField()->getName(), 8 );
    auto currElemChara = _study->getElementaryCharacteristics();
    std::string elemCharaName( " ", 8 );
    if ( currElemChara )
        elemCharaName = std::string( currElemChara->getName(), 0, 8 );

    // Call Fortran WRAPPER
    CALLO_VRCREF( modelName, materialFieldName, elemCharaName, fieldName );

    return field;
}
