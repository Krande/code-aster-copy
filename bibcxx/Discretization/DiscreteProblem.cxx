/**
 * @file DiscreteProblem.cxx
 * @brief Implementation de DiscreteProblem
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

#include "Discretization/DiscreteProblem.h"
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

ElementaryVectorDisplacementRealPtr
DiscreteProblem::computeElementaryDirichletVector( ASTERDOUBLE time ) {
    ElementaryVectorDisplacementRealPtr retour =
        boost::make_shared< ElementaryVectorDisplacementReal >();

    std::string modelName = ljust( _study->getModel()->getName(), 24 );

    JeveuxVectorChar24 jvListOfLoads = _study->getListOfLoads()->getListVector();
    std::string nameLcha = ljust( jvListOfLoads->getName(), 24 );

    JeveuxVectorLong jvInfo = _study->getListOfLoads()->getInformationVector();
    std::string nameInfc = ljust( jvInfo->getName(), 24 );

    std::string typres( "R" );
    std::string resultName( retour->getName() );

    // CORICH appel getres
    CommandSyntax cmdSt( "MECA_STATIQUE" );
    cmdSt.setResult( resultName, "AUCUN" );

    CALLO_VEDIME( modelName, nameLcha, nameInfc, &time, typres, resultName );
    retour->isEmpty( false );

    retour->setListOfLoads( _study->getListOfLoads() );
    return retour;
};

FieldOnNodesRealPtr DiscreteProblem::computeDirichlet( ASTERDOUBLE time ) {
    auto vect_elem = computeElementaryDirichletVector( time );

    return vect_elem->assembleWithLoadFunctions( _study->getDOFNumbering(), time );
};

ElementaryVectorDisplacementRealPtr
DiscreteProblem::computeElementaryDirichletReactionVector( FieldOnNodesRealPtr lagr_curr ) {
    ElementaryVectorDisplacementRealPtr retour =
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

    std::string resultName( retour->getName() );
    const std::string base( "G" );

    std::string lagrName = lagr_curr->getName();

    CALLO_VEBTLA( base, modelName, materName, caraName, lagrName, listLoadsName, resultName );

    retour->isEmpty( false );
    retour->setListOfLoads( listOfLoads );
    retour->build();

    return retour;
};

FieldOnNodesRealPtr DiscreteProblem::computeDirichletReaction( FieldOnNodesRealPtr lagr_curr ) {
    auto vect_elem = computeElementaryDirichletReactionVector( lagr_curr );

    return vect_elem->assemble( _study->getDOFNumbering() );
};

ElementaryVectorDisplacementRealPtr
DiscreteProblem::computeElementaryDualizedDirichletVector( FieldOnNodesRealPtr disp_curr,
                                                           ASTERDOUBLE scaling ) {
    ElementaryVectorDisplacementRealPtr retour =
        boost::make_shared< ElementaryVectorDisplacementReal >();

    std::string modelName = _study->getModel()->getName();
    std::string dispName = disp_curr->getName();

    auto listOfLoads = _study->getListOfLoads();
    std::string listLoadsName = ljust( listOfLoads->getName(), 19 );

    std::string resultName( retour->getName() );
    const std::string base( "G" );
    const ASTERDOUBLE const_scaling = scaling;

    CALLO_VEBUME( modelName, dispName, listLoadsName, resultName, &const_scaling, base );

    retour->isEmpty( false );

    retour->setListOfLoads( listOfLoads );
    return retour;
};

FieldOnNodesRealPtr DiscreteProblem::computeDualizedDirichlet( FieldOnNodesRealPtr disp_curr,
                                                               ASTERDOUBLE scaling ) {
    auto vect_elem = computeElementaryDualizedDirichletVector( disp_curr, scaling );

    auto bume = vect_elem->assemble( _study->getDOFNumbering() );

    if ( _study->getMesh()->isParallel() )
        CALLO_AP_ASSEMBLY_VECTOR( bume->getName() );

    return bume;
};

ElementaryVectorDisplacementRealPtr DiscreteProblem::computeElementaryLaplaceVector() {
    ElementaryVectorDisplacementRealPtr retour =
        boost::make_shared< ElementaryVectorDisplacementReal >();

    ModelPtr curModel = _study->getModel();
    std::string modelName = curModel->getName();
    modelName.resize( 24, ' ' );

    JeveuxVectorChar24 jvListOfLoads = _study->getListOfLoads()->getListVector();
    std::string nameLcha = jvListOfLoads->getName();
    nameLcha.resize( 24, ' ' );

    JeveuxVectorLong jvInfo = _study->getListOfLoads()->getInformationVector();
    std::string nameInfc = jvInfo->getName();
    nameInfc.resize( 24, ' ' );

    std::string blanc( " " );
    const std::string resultName( retour->getName() );

    // CORICH appel getres
    CommandSyntax cmdSt( "MECA_STATIQUE" );
    cmdSt.setResult( resultName, "AUCUN" );

    CALLO_VELAME( modelName, nameLcha, nameInfc, blanc, resultName );
    retour->isEmpty( false );

    retour->setListOfLoads( _study->getListOfLoads() );
    return retour;
};

ElementaryVectorDisplacementRealPtr
DiscreteProblem::computeElementaryNeumannVector( const VectorReal time,
                                                 ExternalStateVariablesBuilderPtr varCom ) {
    if ( time.size() != 3 )
        throw std::runtime_error( "Invalid number of parameter" );

    ElementaryVectorDisplacementRealPtr retour =
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
    if ( varCom != nullptr )
        varCName = varCom->getName() + ".TOUT";
    std::string resultName( retour->getName() );
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
                       retour->getName(), varCName );
    retour->isEmpty( false );

    retour->setListOfLoads( _study->getListOfLoads() );
    return retour;
};

FieldOnNodesRealPtr DiscreteProblem::computeNeumann( const VectorReal time,
                                                     ExternalStateVariablesBuilderPtr varCom ) {
    auto vect_elem = computeElementaryNeumannVector( time, varCom );

    return vect_elem->assembleWithLoadFunctions( _study->getDOFNumbering(), time[0] + time[1] );
};

ElementaryMatrixDisplacementRealPtr
DiscreteProblem::computeElementaryStiffnessMatrix( ASTERDOUBLE time ) {
    ElementaryMatrixDisplacementRealPtr retour =
        boost::make_shared< ElementaryMatrixDisplacementReal >();
    ModelPtr curModel = _study->getModel();
    retour->setModel( curModel );
    MaterialFieldPtr curMater = _study->getMaterialField();
    retour->setMaterialField( curMater );
    auto compor = curMater->getBehaviourField();

    _study->computeListOfLoads();
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
DiscreteProblem::computeElementaryTangentMatrix( ASTERDOUBLE time ) {
    return this->computeElementaryStiffnessMatrix( time );
};

// TODO calcul de la matrice jacobienne pour l'étape de correction de la méthode de Newton
ElementaryMatrixDisplacementRealPtr
DiscreteProblem::computeElementaryJacobianMatrix( ASTERDOUBLE time ) {
    return this->computeElementaryStiffnessMatrix( time );
};

FieldOnNodesRealPtr DiscreteProblem::computeDirichletBC( const ASTERDOUBLE &time ) const {
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

ElementaryVectorDisplacementRealPtr DiscreteProblem::computeElementaryMechanicalLoadsVector() {
    ElementaryVectorDisplacementRealPtr retour( new ElementaryVectorDisplacementReal() );

    CommandSyntax cmdSt( "CALC_VECT_ELEM" );
    cmdSt.setResult( retour->getName(), retour->getType() );

    SyntaxMapContainer dict;
    dict.container["OPTION"] = "CHAR_MECA";
    dict.container["MODELE"] = _study->getModel()->getName();

    if ( _study->getMaterialField() )
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
    cmdSt.define( dict );
    retour->setListOfLoads( listOfLoads );

    try {
        ASTERINTEGER op = 8;
        CALL_EXECOP( &op );
    } catch ( ... ) {
        throw;
    }
    retour->isEmpty( false );

    return retour;
};

SyntaxMapContainer DiscreteProblem::computeMatrixSyntax( const std::string &optionName ) {
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
DiscreteProblem::computeMechanicalMatrix( const std::string &optionName ) {
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

ElementaryMatrixDisplacementRealPtr DiscreteProblem::computeMechanicalDampingMatrix(
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

ElementaryMatrixDisplacementRealPtr DiscreteProblem::computeMechanicalMassMatrix() {
    return computeMechanicalMatrix( "RIGI_MECA" );
};

ElementaryMatrixDisplacementRealPtr DiscreteProblem::computeMechanicalStiffnessMatrix() {
    return computeMechanicalMatrix( "RIGI_MECA" );
};
