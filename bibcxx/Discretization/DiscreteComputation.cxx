/**
 * @file DiscreteComputation.cxx
 * @brief Implementation de DiscreteComputation
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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
#include "Materials/ExternalStateVariablesBuilder.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Modeling/XfemModel.h"
#include "Utilities/Tools.h"

FieldOnNodesRealPtr DiscreteComputation::imposedDisplacement( ASTERDOUBLE currTime ) {

    ElementaryVectorDisplacementRealPtr elemVect =
        boost::make_shared< ElementaryVectorDisplacementReal >();

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
    elemVect->build();
    elemVect->isEmpty( false );
    elemVect->setListOfLoads( listOfLoads );

    // Assemble
    return elemVect->assembleWithLoadFunctions( _study->getDOFNumbering(), currTime );
};

FieldOnNodesRealPtr DiscreteComputation::dualReaction( FieldOnNodesRealPtr lagr_curr ) {

    ElementaryVectorDisplacementRealPtr elemVect =
        boost::make_shared< ElementaryVectorDisplacementReal >();

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
    elemVect->build();
    elemVect->isEmpty( false );
    elemVect->setListOfLoads( listOfLoads );

    // Assemble
    return elemVect->assemble( _study->getDOFNumbering() );
};

FieldOnNodesRealPtr DiscreteComputation::dualDisplacement( FieldOnNodesRealPtr disp_curr,
                                                           ASTERDOUBLE scaling ) {

    ElementaryVectorDisplacementRealPtr elemVect =
        boost::make_shared< ElementaryVectorDisplacementReal >();

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
    elemVect->build();
    elemVect->isEmpty( false );
    elemVect->setListOfLoads( listOfLoads );

    // Assemble
    FieldOnNodesRealPtr bume = elemVect->assemble( _study->getDOFNumbering() );

    if ( _study->getMesh()->isParallel() )
        CALLO_AP_ASSEMBLY_VECTOR( bume->getName() );

    return bume;
};

FieldOnNodesRealPtr DiscreteComputation::neumann( const VectorReal timeParameters ) {

    ElementaryVectorDisplacementRealPtr elemVect =
        boost::make_shared< ElementaryVectorDisplacementReal >();

    if ( _study->getModel()->isThermal() ) {
        AS_ABORT( "Not implemented for thermic" );
    };

    if ( timeParameters.size() != 3 )
        throw std::runtime_error( "Invalid number of parameter" );
    const ASTERDOUBLE &currTime = timeParameters[0];

    // Prepare loads
    auto listOfLoads = _study->getListOfLoads();
    JeveuxVectorChar24 listOfLoadsList = listOfLoads->getListVector();
    JeveuxVectorLong listOfLoadsInfo = listOfLoads->getInformationVector();
    std::string nameLcha = ljust( listOfLoadsList->getName(), 24 );
    std::string nameInfc = ljust( listOfLoadsInfo->getName(), 24 );

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( _study->getModel()->getName(), 24 );
    std::string materName = ljust( _study->getMaterialField()->getName(), 24 );
    std::string codmaName =
        ljust( _study->getCodedMaterial()->getCodedMaterialField()->getName(), 24 );
    auto currElemChara = _study->getElementaryCharacteristics();
    std::string caraName( " " );
    if ( currElemChara )
        caraName = currElemChara->getName();
    caraName.resize( 24, ' ' );
    std::string varcName( " " );
    auto varCom = _study->getExternalStateVariables();
    if ( varCom != nullptr ) {
        varCom->build( currTime );
        varcName = varCom->getName() + ".TOUT";
    }
    varcName.resize( 24, ' ' );
    std::string vectElemName = ljust( elemVect->getName(), 24 );
    std::string stop( "S" );

    // Wrapper FORTRAN
    CALLO_VECHME_WRAP( stop, modelName, nameLcha, nameInfc, &currTime, caraName, materName,
                       codmaName, vectElemName, varcName );

    // Construct vect_elem object
    elemVect->build();
    elemVect->isEmpty( false );
    elemVect->setListOfLoads( _study->getListOfLoads() );

    // Assemble
    return elemVect->assembleWithLoadFunctions( _study->getDOFNumbering(),
                                                timeParameters[0] + timeParameters[1] );
};

FieldOnNodesRealPtr DiscreteComputation::dirichletBC( const ASTERDOUBLE &time ) const {

    auto dofNume = _study->getDOFNumbering();
    FieldOnNodesRealPtr vectAsse = boost::make_shared< FieldOnNodesReal >( dofNume );

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

        return boost::make_shared< FieldOnNodesReal >( diri_impo );
    }

    // Construct vect_asse object
    FieldOnNodesRealPtr vectAsse = boost::make_shared< FieldOnNodesReal >( dofNume );
    vectAsse->setValues( 0.0 );
    vectAsse->build();

    return vectAsse;
};

FieldOnNodesRealPtr DiscreteComputation::externalStateVariables( const ASTERDOUBLE &time ) {
    auto varCom = _study->getExternalStateVariables();

    varCom->build( time );

    return varCom->computeExternalStateVariablesLoad( _study->getDOFNumbering() );
};

ElementaryMatrixDisplacementRealPtr
DiscreteComputation::elasticStiffnessMatrix( ASTERDOUBLE time ) {

    ElementaryMatrixDisplacementRealPtr elemMatr =
        boost::make_shared< ElementaryMatrixDisplacementReal >();

    // Get main parameters
    const std::string option( "RIGI_MECA" );
    int nh = 0;
    ModelPtr currModel = _study->getModel();
    MaterialFieldPtr currMater = _study->getMaterialField();
    CodedMaterialPtr currCodedMater = _study->getCodedMaterial();
    ElementaryCharacteristicsPtr currElemChara = _study->getElementaryCharacteristics();
    ExternalStateVariablesBuilderPtr currExteVari = _study->getExternalStateVariables();

    // Compute external state variables
    if ( currExteVari->hasExternalStateVariables() ||
         currMater->hasExternalStateVariables( "NEUT1" ) ||
         currMater->hasExternalStateVariables( "GEOM" ) ) {
        currExteVari->build( time );
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
    _calcul->setModel( currModel );
    elemMatr->prepareCompute( option );

    // Add input fields
    _calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    _calcul->addInputField( "PVARCPR", currExteVari->getExternalStateVariablesField() );
    _calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
    _calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );
    if ( currElemChara ) {
        _calcul->addElementaryCharacteristicsField( currElemChara );
    }
    _calcul->addFourierModeField( nh );
    _calcul->addTimeField( time );
    if ( currModel->existsXfem() ) {
        XfemModelPtr currXfemModel = currModel->getXfemModel();
        _calcul->addXFEMField( currXfemModel );
    }

    // Add output elementary terms
    ElementaryTermRealPtr resuElemSyme =
        boost::make_shared< ElementaryTermReal >( elemMatr->generateNameOfElementaryTerm() );
    _calcul->addOutputElementaryTerm( "PMATUUR", resuElemSyme );
    ElementaryTermRealPtr resuElemNSym =
        boost::make_shared< ElementaryTermReal >( elemMatr->generateNameOfElementaryTerm() );
    _calcul->addOutputElementaryTerm( "PMATUNS", resuElemNSym );

    // Compute elementary matrices for rigidity
    if ( currModel->existsFiniteElement() ) {
        _calcul->compute();
        if ( _calcul->hasOutputElementaryTerm( "PMATUUR" ) )
            elemMatr->addElementaryTerm( _calcul->getOutputElementaryTerm( "PMATUUR" ) );
        if ( _calcul->hasOutputElementaryTerm( "PMATUNS" ) )
            elemMatr->addElementaryTerm( _calcul->getOutputElementaryTerm( "PMATUNS" ) );
    };

    // Compute elementary matrices for dual BC
    DiscreteComputation::baseComputeMechanicalDualBCMatrix( _calcul, elemMatr );

    elemMatr->isEmpty( false );
    return elemMatr;
};

SyntaxMapContainer DiscreteComputation::computeMatrixSyntax( const std::string &option ) {
    SyntaxMapContainer dict;

    // Definition du mot cle simple MODELE
    if ( ( !_study->getModel() ) || _study->getModel()->isEmpty() )
        throw std::runtime_error( "Model is empty" );
    dict.container["MODELE"] = _study->getModel()->getName();

    // Definition du mot cle simple CHAM_MATER
    if ( !_study->getMaterialField() )
        throw std::runtime_error( "Material is empty" );
    dict.container["CHAM_MATER"] = _study->getMaterialField()->getName();

    auto listOfLoads = _study->getListOfLoads();

    const auto listOfMechanicalLoadReal = listOfLoads->getMechanicalLoadsReal();
    if ( listOfMechanicalLoadReal.size() != 0 ) {
        VectorString tmp;
        for ( const auto curIter : listOfMechanicalLoadReal )
            tmp.push_back( curIter->getName() );
        dict.container["CHARGE"] = tmp;
    }

    const auto listOfMechanicalLoadFunction = listOfLoads->getMechanicalLoadsFunction();
    if ( listOfMechanicalLoadFunction.size() != 0 ) {
        VectorString tmp;
        for ( const auto curIter : listOfMechanicalLoadFunction )
            tmp.push_back( curIter->getName() );
        dict.container["CHARGE"] = tmp;
    }
#ifdef ASTER_HAVE_MPI
    auto listParaMecaLoadReal = listOfLoads->getParallelMechanicalLoadsReal();
    if ( listParaMecaLoadReal.size() != 0 ) {
        VectorString tmp;
        for ( const auto curIter : listParaMecaLoadReal )
            tmp.push_back( curIter->getName() );
        dict.container["CHARGE"] = tmp;
    }

    auto listParaMecaLoadFunction = listOfLoads->getParallelMechanicalLoadsFunction();
    if ( listParaMecaLoadFunction.size() != 0 ) {
        VectorString tmp;
        for ( const auto curIter : listParaMecaLoadFunction )
            tmp.push_back( curIter->getName() );
        dict.container["CHARGE"] = tmp;
    }
#endif /* ASTER_HAVE_MPI */

    // Definition du mot cle simple OPTION
    dict.container["OPTION"] = option;

    return dict;
};

ElementaryMatrixDisplacementRealPtr DiscreteComputation::computeMechanicalDampingMatrix(
    const ElementaryMatrixDisplacementRealPtr &rigiMatrElem,
    const ElementaryMatrixDisplacementRealPtr &massMatrElem ) {
    ElementaryMatrixDisplacementRealPtr elemMatr( new ElementaryMatrixDisplacementReal() );
    elemMatr->setModel( rigiMatrElem->getModel() );

    // Definition du bout de fichier de commande correspondant a CALC_MATR_ELEM
    CommandSyntax cmdSt( "CALC_MATR_ELEM" );
    cmdSt.setResult( elemMatr->getName(), elemMatr->getType() );

    SyntaxMapContainer dict = computeMatrixSyntax( "AMOR_MECA" );
    dict.container["RIGI_MECA"] = rigiMatrElem->getName();
    dict.container["MASS_MECA"] = massMatrElem->getName();

    cmdSt.define( dict );
    try {
        ASTERINTEGER op = 9;
        CALL_EXECOP( &op );
    } catch ( ... ) {
        throw;
    }
    elemMatr->isEmpty( false );
    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixDisplacementRealPtr DiscreteComputation::computeMechanicalStiffnessMatrix() {

    ASTERDOUBLE time = 0.0;
    return elasticStiffnessMatrix( time );
};

ElementaryMatrixDisplacementRealPtr DiscreteComputation::computeMechanicalDualBCMatrix() {

    ElementaryMatrixDisplacementRealPtr elemMatr =
        boost::make_shared< ElementaryMatrixDisplacementReal >();

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
    CalculPtr _calcul = boost::make_shared< Calcul >( option );
    elemMatr->prepareCompute( option );

    // Compute elementary matrices
    DiscreteComputation::baseComputeMechanicalDualBCMatrix( _calcul, elemMatr );

    elemMatr->isEmpty( false );
    return elemMatr;
};

void DiscreteComputation::baseComputeMechanicalDualBCMatrix(
    CalculPtr &calcul, ElementaryMatrixDisplacementRealPtr &elemMatr ) {

    // Prepare loads
    const auto &_listOfLoads = _study->getListOfLoads();
    ListMecaLoadReal mecaLoadReal = _listOfLoads->getMechanicalLoadsReal();

    // Select option
    calcul->setOption( "MECA_DDLM_R" );
    calcul->clearInputs();
    calcul->clearOutputs();

    for ( const auto &curIter : mecaLoadReal ) {
        FiniteElementDescriptorPtr FEDesc = curIter->getFiniteElementDescriptor();
        ConstantFieldOnCellsRealPtr field =
            curIter->getMechanicalLoadDescription()->getMultiplicativeField();
        if ( field && field->exists() && FEDesc ) {
            std::string elemTermName = elemMatr->generateNameOfElementaryTerm();
            ElementaryTermRealPtr resuElem =
                boost::make_shared< ElementaryTermReal >( elemTermName );
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( FEDesc );
            calcul->addInputField( "PDDLMUR", field );
            calcul->addOutputElementaryTerm( "PMATUUR", resuElem );
            calcul->compute();
            elemMatr->addElementaryTerm( resuElem );
        }
    }

    ListMecaLoadFunction mecaLoadFunc = _listOfLoads->getMechanicalLoadsFunction();
    for ( const auto &curIter : mecaLoadFunc ) {
        FiniteElementDescriptorPtr FEDesc = curIter->getFiniteElementDescriptor();
        ConstantFieldOnCellsRealPtr field =
            curIter->getMechanicalLoadDescription()->getMultiplicativeField();
        if ( field && field->exists() && FEDesc ) {
            std::string elemTermName = elemMatr->generateNameOfElementaryTerm();
            ElementaryTermRealPtr resuElem =
                boost::make_shared< ElementaryTermReal >( elemTermName );
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( FEDesc );
            calcul->addInputField( "PDDLMUR", field );
            calcul->addOutputElementaryTerm( "PMATUUR", resuElem );
            calcul->compute();
            elemMatr->addElementaryTerm( resuElem );
        }
    }
};
