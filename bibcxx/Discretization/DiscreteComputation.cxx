/**
 * @file DiscreteComputation.cxx
 * @brief Implementation de DiscreteComputation
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

#include <iostream>
#include <string>

#include "Discretization/DiscreteComputation.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Materials/ExternalStateVariablesBuilder.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Numbering/ParallelDOFNumbering.h"
#include "Utilities/Tools.h"
#include "aster_fort_calcul.h"
#include "aster_fort_superv.h"

/* person_in_charge: nicolas.sellenet at edf.fr */

FieldOnNodesRealPtr DiscreteComputation::imposedDisplacement( ASTERDOUBLE time ) {
    ElementaryVectorDisplacementRealPtr vect_elem =
        boost::make_shared< ElementaryVectorDisplacementReal >();

    std::string modelName = ljust( _study->getModel()->getName(), 24 );

    JeveuxVectorChar24 jvListOfLoads = _study->getListOfLoads()->getListVector();
    std::string nameLcha = ljust( jvListOfLoads->getName(), 24 );

    JeveuxVectorLong jvInfo = _study->getListOfLoads()->getInformationVector();
    std::string nameInfc = ljust( jvInfo->getName(), 24 );

    std::string typres( "R" );
    std::string resultName( vect_elem->getName() );

    // CORICH appel getres
    CommandSyntax cmdSt( "MECA_STATIQUE" );
    cmdSt.setResult( resultName, "AUCUN" );

    CALLO_VEDIME( modelName, nameLcha, nameInfc, &time, typres, resultName );

    vect_elem->isEmpty( false );
    vect_elem->setListOfLoads( _study->getListOfLoads() );
    cmdSt.free();

    return vect_elem->assembleWithLoadFunctions( _study->getDOFNumbering(), time );
};

FieldOnNodesRealPtr DiscreteComputation::dualReaction( FieldOnNodesRealPtr lagr_curr ) {

    ElementaryVectorDisplacementRealPtr vect_elem =
        boost::make_shared< ElementaryVectorDisplacementReal >();

    std::string modelName = ljust( _study->getModel()->getName(), 24 );
    std::string materName = ljust( _study->getMaterialField()->getName(), 24 );

    auto curCaraElem = _study->getElementaryCharacteristics();
    std::string caraName( " " );
    if ( curCaraElem )
        caraName = curCaraElem->getName();
    caraName.resize( 24, ' ' );

    auto listOfLoads = _study->getListOfLoads();
    std::string listLoadsName = ljust( listOfLoads->getName(), 19 );

    std::string resultName( vect_elem->getName() );
    const std::string base( "G" );

    std::string lagrName = lagr_curr->getName();

    CALLO_VEBTLA( base, modelName, materName, caraName, lagrName, listLoadsName, resultName );

    vect_elem->isEmpty( false );
    vect_elem->setListOfLoads( listOfLoads );
    vect_elem->build();

    return vect_elem->assemble( _study->getDOFNumbering() );
};

FieldOnNodesRealPtr DiscreteComputation::dualDisplacement( FieldOnNodesRealPtr disp_curr,
                                                           ASTERDOUBLE scaling ) {

    ElementaryVectorDisplacementRealPtr vect_elem =
        boost::make_shared< ElementaryVectorDisplacementReal >();

    std::string modelName = _study->getModel()->getName();
    std::string dispName = disp_curr->getName();

    auto listOfLoads = _study->getListOfLoads();
    std::string listLoadsName = ljust( listOfLoads->getName(), 19 );

    std::string resultName( vect_elem->getName() );
    const std::string base( "G" );
    const ASTERDOUBLE const_scaling = scaling;

    CALLO_VEBUME( modelName, dispName, listLoadsName, resultName, &const_scaling, base );

    vect_elem->isEmpty( false );

    vect_elem->setListOfLoads( listOfLoads );

    FieldOnNodesRealPtr bume = vect_elem->assemble( _study->getDOFNumbering() );

    if ( _study->getMesh()->isParallel() )
        CALLO_AP_ASSEMBLY_VECTOR( bume->getName() );

    return bume;
};

FieldOnNodesRealPtr DiscreteComputation::neumann( const VectorReal time ) {

    if ( time.size() != 3 )
        throw std::runtime_error( "Invalid number of parameter" );

    ElementaryVectorDisplacementRealPtr vect_elem =
        boost::make_shared< ElementaryVectorDisplacementReal >();
    const auto &curCodedMater = _study->getCodedMaterial()->getCodedMaterialField();
    const auto &curMater = _study->getCodedMaterial()->getMaterialField();

    ModelPtr curModel = _study->getModel();
    std::string modelName = curModel->getName();

    JeveuxVectorChar24 jvListOfLoads = _study->getListOfLoads()->getListVector();
    std::string nameLcha = jvListOfLoads->getName();

    JeveuxVectorLong jvInfo = _study->getListOfLoads()->getInformationVector();
    std::string nameInfc = jvInfo->getName();

    const ASTERDOUBLE &inst = time[0];
    std::string stop( "S" );
    std::string blanc( "        " );
    std::string varCName( blanc );
    auto varCom = _study->getExternalStateVariables();
    if ( varCom != nullptr )
    {
        varCom->build(inst);
        varCName = varCom->getName() + ".TOUT";
    }
    std::string resultName( vect_elem->getName() );
    std::string materName( curMater->getName() + "                " );
    std::string codmaName( curCodedMater->getName() + "                " );

    std::string caraName( blanc );
    const auto &caraElem = _study->getElementaryCharacteristics();
    if ( caraElem != nullptr )
        caraName = caraElem->getName();

    // CORICH appel getres
    CommandSyntax cmdSt( "MECA_STATIQUE" );
    cmdSt.setResult( resultName, "AUCUN" );

    CALLO_VECHME_WRAP( stop, modelName, nameLcha, nameInfc, &inst, caraName, materName, codmaName,
                       vect_elem->getName(), varCName );

    vect_elem->isEmpty( false );
    vect_elem->setListOfLoads( _study->getListOfLoads() );
    vect_elem->build();
    cmdSt.free();

    return vect_elem->assembleWithLoadFunctions( _study->getDOFNumbering(), time[0] + time[1] );
};

FieldOnNodesRealPtr DiscreteComputation::dirichletBC( const ASTERDOUBLE &time ) const {
    const auto &_listOfLoad = _study->getListOfLoads();
    const auto &list = _listOfLoad->getListVector();
    const auto &loadInformations = _listOfLoad->getInformationVector();
    const auto &listOfFunctions = _listOfLoad->getListOfFunctions();
    if ( _listOfLoad->isEmpty() )
        _listOfLoad->build( _study->getModel() );
    //         throw std::runtime_error( "ListOfLoads is empty" );

    FieldOnNodesRealPtr retour = boost::make_shared< FieldOnNodesReal >();
    std::string resuName = retour->getName();
    std::string dofNumName = _study->getDOFNumbering()->getName();

    std::string lLoadName = list->getName();
    lLoadName.resize( 24, ' ' );
    std::string infLoadName = loadInformations->getName();
    infLoadName.resize( 24, ' ' );
    std::string funcLoadName = listOfFunctions->getName();
    funcLoadName.resize( 24, ' ' );
    std::string base( "G" );

    CALLO_ASCAVC_WRAP( lLoadName, infLoadName, funcLoadName, dofNumName, &time, resuName, base );

    retour->setDOFNumbering( _study->getDOFNumbering() );
    retour->build();

    return retour;
};

FieldOnNodesRealPtr
DiscreteComputation::incrementalDirichletBC( const ASTERDOUBLE &time,
                                             const FieldOnNodesRealPtr disp_curr ) const {
    auto dofNume = _study->getDOFNumbering();

    if ( dofNume->hasDirichletBC() ) {
        auto diri_curr = dirichletBC( time );
        auto diri_impo = *(diri_curr) - *(disp_curr);

        // Set to zero terms not imposed
        auto eliminatedDofs = dofNume->getDirichletBCDOFs();
        auto nbElimination = eliminatedDofs.size();

        for ( ASTERINTEGER ieq = 0; ieq < nbElimination; ieq++ ) {
            if ( eliminatedDofs[ieq] == 0 )
                diri_impo[ieq] = 0.0;
        }

        return boost::make_shared< FieldOnNodesReal >(diri_impo);
    }

    FieldOnNodesRealPtr diri_impo = boost::make_shared< FieldOnNodesReal >( dofNume );
    diri_impo->setValues( 0.0 );
    diri_impo->build();

    return diri_impo;
};

FieldOnNodesRealPtr
DiscreteComputation::externalStateVariables( const ASTERDOUBLE &time )  {
    auto varCom = _study->getExternalStateVariables();

    varCom->build(time);

    return varCom->computeExternalStateVariablesLoad(_study->getDOFNumbering());
};

ElementaryMatrixDisplacementRealPtr
DiscreteComputation::elasticStiffnessMatrix( ASTERDOUBLE time ) {
    ElementaryMatrixDisplacementRealPtr retour =
        boost::make_shared< ElementaryMatrixDisplacementReal >();
    ModelPtr curModel = _study->getModel();
    retour->setModel( curModel );
    MaterialFieldPtr curMater = _study->getMaterialField();
    retour->setMaterialField( curMater );
    auto compor = curMater->getBehaviourField();

    JeveuxVectorChar24 jvListOfLoads = _study->getListOfLoads()->getListVector();
    jvListOfLoads->updateValuePointer();
    ASTERINTEGER nbLoad = jvListOfLoads->size();

    std::string blanc( 24, ' ' );
    std::string modelName = curModel->getName();
    modelName.resize( 24, ' ' );

    const auto &codedMater = _study->getCodedMaterial()->getCodedMaterialField();

    std::string caraName( blanc );
    const auto &caraElem = _study->getElementaryCharacteristics();
    if ( caraElem != nullptr )
        caraName = caraElem->getName();

    std::string materName = curMater->getName();
    materName.resize( 24, ' ' );
    std::string coMatName = codedMater->getName();
    coMatName.resize( 24, ' ' );

    // MERIME appel getres
    CommandSyntax cmdSt( "MECA_STATIQUE" );
    cmdSt.setResult( "AUCUN", "AUCUN" );
    SyntaxMapContainer dict;
    dict.container["INFO"] = (ASTERINTEGER)1;
    cmdSt.define( dict );

    ASTERINTEGER nh = 0;

    CALLO_MERIME_WRAP( modelName, &nbLoad, *( jvListOfLoads->getDataPtr() ), materName, coMatName,
                       caraName, &time, compor->getName(), retour->getName(), &nh,
                       JeveuxMemoryTypesNames[0] );

    retour->build();
    retour->isEmpty( false );
    return retour;
};

// TODO calcul de la matrice tangente pour l'étape de prédiction de la méthode de Newton
ElementaryMatrixDisplacementRealPtr
DiscreteComputation::computeElementaryTangentMatrix( ASTERDOUBLE time ) {
    return this->elasticStiffnessMatrix( time );
};

// TODO calcul de la matrice jacobienne pour l'étape de correction de la méthode de Newton
ElementaryMatrixDisplacementRealPtr
DiscreteComputation::computeElementaryJacobianMatrix( ASTERDOUBLE time ) {
    return this->elasticStiffnessMatrix( time );
};

SyntaxMapContainer DiscreteComputation::computeMatrixSyntax( const std::string &optionName ) {
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
    dict.container["OPTION"] = optionName;

    return dict;
};

ElementaryMatrixDisplacementRealPtr
DiscreteComputation::computeMechanicalMatrix( const std::string &optionName ) {
    ElementaryMatrixDisplacementRealPtr retour( new ElementaryMatrixDisplacementReal() );
    retour->setModel( _study->getModel() );

    // Definition du bout de fichier de commande correspondant a CALC_MATR_ELEM
    CommandSyntax cmdSt( "CALC_MATR_ELEM" );
    cmdSt.setResult( retour->getName(), retour->getType() );

    SyntaxMapContainer dict = computeMatrixSyntax( optionName );

    cmdSt.define( dict );
    try {
        ASTERINTEGER op = 9;
        CALL_EXECOP( &op );
    } catch ( ... ) {
        throw;
    }
    retour->build();
    retour->isEmpty( false );

    return retour;
};

ElementaryMatrixDisplacementRealPtr DiscreteComputation::computeMechanicalDampingMatrix(
    const ElementaryMatrixDisplacementRealPtr &rigidity,
    const ElementaryMatrixDisplacementRealPtr &mass ) {
    ElementaryMatrixDisplacementRealPtr retour( new ElementaryMatrixDisplacementReal() );
    retour->setModel( rigidity->getModel() );

    // Definition du bout de fichier de commande correspondant a CALC_MATR_ELEM
    CommandSyntax cmdSt( "CALC_MATR_ELEM" );
    cmdSt.setResult( retour->getName(), retour->getType() );

    SyntaxMapContainer dict = computeMatrixSyntax( "AMOR_MECA" );
    dict.container["RIGI_MECA"] = rigidity->getName();
    dict.container["MASS_MECA"] = mass->getName();

    cmdSt.define( dict );
    try {
        ASTERINTEGER op = 9;
        CALL_EXECOP( &op );
    } catch ( ... ) {
        throw;
    }
    retour->isEmpty( false );
    retour->build();
    return retour;
};

ElementaryMatrixDisplacementRealPtr DiscreteComputation::computeMechanicalMassMatrix() {
    return computeMechanicalMatrix( "RIGI_MECA" );
};

ElementaryMatrixDisplacementRealPtr DiscreteComputation::computeMechanicalStiffnessMatrix() {
    return computeMechanicalMatrix( "RIGI_MECA" );
};
