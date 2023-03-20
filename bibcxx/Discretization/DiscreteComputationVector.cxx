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
#include "Supervis/CommandSyntax.h"
#include "Utilities/Tools.h"

FieldOnNodesRealPtr DiscreteComputation::getImposedDualBC( const ASTERDOUBLE time,
                                                           const ASTERDOUBLE time_delta,
                                                           const ASTERDOUBLE time_theta ) const {

    bool has_load = false;

    auto elemVect = std::make_shared< ElementaryVectorReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

    if ( _phys_problem->getModel()->isThermal() ) {
        has_load = this->addTherImposedTerms( elemVect, time, time_delta, time_theta );
    } else if ( _phys_problem->getModel()->isMechanical() ) {
        has_load = this->addMecaImposedTerms( elemVect, time );
    } else {
        AS_ASSERT( false );
    };

    if ( has_load ) {
        auto FEDs = _phys_problem->getListOfLoads()->getFiniteElementDescriptors();
        FEDs.push_back( _phys_problem->getModel()->getFiniteElementDescriptor() );
        elemVect->build( FEDs );
        return elemVect->assembleWithLoadFunctions( _phys_problem->getDOFNumbering(), time );
    } else {
        FieldOnNodesRealPtr vectAsse =
            std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
        vectAsse->setValues( 0.0 );
        return vectAsse;
    }
};

FieldOnNodesRealPtr
DiscreteComputation::getNeumannForces( const ASTERDOUBLE time, const ASTERDOUBLE time_delta,
                                       const ASTERDOUBLE time_theta,
                                       const FieldOnNodesRealPtr _previousPrimalField ) const {

    bool has_load = false;

    auto elemVect = std::make_shared< ElementaryVectorReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

    if ( _phys_problem->getModel()->isThermal() ) {
        has_load = this->addTherNeumannTerms( elemVect, time, time_delta, time_theta,
                                              _previousPrimalField );
    } else if ( _phys_problem->getModel()->isMechanical() ) {
        has_load = this->addMecaNeumannTerms( elemVect, time, time_delta, time_theta );
    } else {
        AS_ASSERT( false );
    };
    if ( has_load ) {
        auto FEDs = _phys_problem->getListOfLoads()->getFiniteElementDescriptors();
        FEDs.push_back( _phys_problem->getModel()->getFiniteElementDescriptor() );
        elemVect->build( FEDs );
        return elemVect->assembleWithLoadFunctions( _phys_problem->getDOFNumbering(), time );
    } else {
        FieldOnNodesRealPtr vectAsse =
            std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
        vectAsse->setValues( 0.0 );
        vectAsse->build();
        return vectAsse;
    }
};

FieldOnNodesRealPtr DiscreteComputation::getDualForces( FieldOnNodesRealPtr lagr_curr ) const {

    auto elemVect = std::make_shared< ElementaryVectorReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

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
    auto FEDs = _phys_problem->getListOfLoads()->getFiniteElementDescriptors();
    FEDs.push_back( _phys_problem->getModel()->getFiniteElementDescriptor() );
    elemVect->build( FEDs );

    // Assemble
    return elemVect->assemble( _phys_problem->getDOFNumbering() );
};

FieldOnNodesRealPtr DiscreteComputation::getDualDisplacement( FieldOnNodesRealPtr disp_curr,
                                                              ASTERDOUBLE scaling ) const {

    auto elemVect = std::make_shared< ElementaryVectorReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

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
    auto FEDs = _phys_problem->getListOfLoads()->getFiniteElementDescriptors();
    FEDs.push_back( _phys_problem->getModel()->getFiniteElementDescriptor() );
    elemVect->build( FEDs );

    // Assemble
    FieldOnNodesRealPtr bume = elemVect->assemble( _phys_problem->getDOFNumbering() );

    if ( _phys_problem->getMesh()->isParallel() )
        CALLO_AP_ASSEMBLY_VECTOR( bume->getName() );

    return bume;
};

template < typename T >
std::shared_ptr< FieldOnNodes< T > >
DiscreteComputation::_getDirichletBC( const ASTERDOUBLE time ) const {

    auto dofNume = _phys_problem->getDOFNumbering();
    auto vectAsse = std::make_shared< FieldOnNodes< T > >( dofNume );

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
    std::string dofNumName = dofNume->getName();
    std::string base( "G" );

    // Wrapper FORTRAN
    CALLO_ASCAVC_WRAP( nameLcha, nameInfc, nameFcha, dofNumName, &time, vectAsseName, base );

    // Construct vect_asse object
    vectAsse->build();

    // Assemble
    return vectAsse;
};

FieldOnNodesRealPtr DiscreteComputation::getMechanicalDirichletBC( const ASTERDOUBLE time ) const {
    return this->_getDirichletBC< ASTERDOUBLE >( time );
}

FieldOnNodesRealPtr DiscreteComputation::getThermalDirichletBC( const ASTERDOUBLE time ) const {
    return this->_getDirichletBC< ASTERDOUBLE >( time );
}

FieldOnNodesComplexPtr DiscreteComputation::getAcousticDirichletBC( const ASTERDOUBLE time ) const {
    CommandSyntax cmdSt( "CALC_CHAR_CINE" );
    SyntaxMapContainer dict;
    ListSyntaxMapContainer listeExcit;
    listeExcit.push_back( dict );
    SyntaxMapContainer dict2;

    dict.container["EXCIT"] = listeExcit;
    cmdSt.define( dict );

    return this->_getDirichletBC< ASTERCOMPLEX >( time );
}

FieldOnNodesRealPtr
DiscreteComputation::getIncrementalDirichletBC( const ASTERDOUBLE &time,
                                                const FieldOnNodesRealPtr disp_curr ) const {
    auto dofNume = _phys_problem->getDOFNumbering();

    const auto listOfLoads = _phys_problem->getListOfLoads();
    if ( listOfLoads->hasDirichletBC() ) {
        auto diri_curr = this->_getDirichletBC< ASTERDOUBLE >( time );
        auto diri_impo = *( diri_curr ) - *( disp_curr );
        diri_impo.updateValuePointers();

        // Set to zero terms not imposed
        auto eliminatedDofs = _phys_problem->getDirichletBCDOFs();
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

    return vectAsse;
};

FieldOnNodesRealPtr
DiscreteComputation::getExternalStateVariablesForces( const ASTERDOUBLE time ) const {

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currExternVarRefe = _phys_problem->getReferenceExternalStateVariables();

    // Some checks
    AS_ASSERT( currMater );
    AS_ASSERT( currMater->hasExternalStateVariableForLoad() );
    if ( currMater->hasExternalStateVariableWithReference() ) {
        AS_ASSERT( currExternVarRefe );
    }

    // Main object
    CalculPtr calcul = std::make_unique< Calcul >( "CHAR_VARC" );

    // Create specific output field for XFEM
    FieldOnCellsRealPtr sigmXfem;
    if ( currModel->existsXfem() ) {
        const std::string option = "SIEF_ELGA";
        const std::string paraName = "PCONTRR";
        sigmXfem = std::make_shared< FieldOnCellsReal >( currModel, option, paraName );
    }

    // Create elementary vectors
    auto elemVect = std::make_shared< ElementaryVectorReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );
    elemVect->prepareCompute( "CHAR_VARC" );

    int nbExternVar = static_cast< int >( externVarEnumInt::NumberOfExternVarTypes );
    for ( auto iExternVar = 0; iExternVar < nbExternVar; iExternVar++ ) {
        externVarEnumInt numeExternVar = static_cast< externVarEnumInt >( iExternVar );
        if ( currMater->hasExternalStateVariable( numeExternVar ) &&
             ExternalVariableTraits::externVarHasStrain( numeExternVar ) ) {
            const auto option = ExternalVariableTraits::getExternVarOption( numeExternVar );
            calcul->setOption( option );
            calcul->setModel( currModel );
            calcul->clearInputs();
            calcul->clearOutputs();

            // Add input fields
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
            calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );
            if ( currElemChara ) {
                calcul->addElementaryCharacteristicsField( currElemChara );
            }
            calcul->addInputField( "PVARCPR", _phys_problem->getExternalStateVariables( time ) );
            calcul->addTimeField( "PTEMPSR", time );
            if ( currExternVarRefe ) {
                calcul->addInputField( "PVARCRR", currExternVarRefe );
            }
            if ( currModel->existsXfem() ) {
                XfemModelPtr currXfemModel = currModel->getXfemModel();
                calcul->addXFEMField( currXfemModel );
                calcul->addOutputField( "PCONTRT", sigmXfem );
            }
            calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ) );
            }
        }
    }
    // Build elementary vectors
    elemVect->build();

    // Assemble
    return elemVect->assemble( _phys_problem->getDOFNumbering() );
}
