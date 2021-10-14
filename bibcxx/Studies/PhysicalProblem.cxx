/**
 * @file PhysicalProblem.cxx
 * @brief Implementation de PhysicalProblem
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

#include "Studies/PhysicalProblem.h"
#include "Supervis/Exceptions.h"

PhysicalProblem::PhysicalProblem( const ModelPtr curModel, const MaterialFieldPtr curMat,
                                  const ElementaryCharacteristicsPtr cara )
    : _model( curModel ), _mesh( curModel->getMesh() ), _materialField( curMat ),
      _listOfLoads( boost::make_shared< ListOfLoads >(_model) ), _elemChara( cara ),
      _codedMater( boost::make_shared< CodedMaterial >( _materialField, _model ) ),
      _varCom( boost::make_shared< ExternalStateVariablesBuilder >( _model, _materialField,
                                                                    _elemChara, _codedMater ) ),
      _behavProp( boost::make_shared< BehaviourProperty >( _model, _materialField ) ) {
    // Add checks
    if ( _elemChara ) {
        if ( _model->getName() != _elemChara->getModel()->getName() )
            raiseAsterError( "Inconsistent model" );
    }

    if ( _mesh->getName() != _materialField->getMesh()->getName() )
        raiseAsterError( "Inconsistent mesh" );

// create dofNume
#ifdef ASTER_HAVE_MPI
    if ( _mesh->isParallel() )
        _dofNume = boost::make_shared< ParallelDOFNumbering >();
    else
#endif /* ASTER_HAVE_MPI */
        _dofNume = boost::make_shared< DOFNumbering >();

    _dofNume->setModel( _model );
    _dofNume->setListOfLoads( _listOfLoads );
};

void PhysicalProblem::setDOFNumbering( const BaseDOFNumberingPtr dofNume ) {
    if ( dofNume->getMesh()->getName() != _mesh->getName() )
        raiseAsterError( "Incompatible meshes" );

    auto model = dofNume->getModel();
    if ( model && model->getName() != _model->getName() )
        raiseAsterError( "Incompatible models" );

    auto listOfLoads = dofNume->getListOfLoads();
    if ( listOfLoads && listOfLoads->getName() != _listOfLoads->getName() )
        raiseAsterError( "Incompatible list of loads" );

    _dofNume = dofNume;
};

bool PhysicalProblem::computeDOFNumbering() { return _dofNume->computeNumbering(); };

void PhysicalProblem::computeBehaviourProperty( PyObject *keywords, const std::string &initialState,
                                                const std::string &implex,
                                                const ASTERINTEGER verbosity ) {
    // Create object for behaviour
    _behavProp->setInitialState( initialState == "OUI" );
    _behavProp->setImplex( implex == "OUI" );
    _behavProp->setVerbosity( verbosity > 1 );

    // Check input PyObject
    if ( !PyDict_Check( keywords ) && !PyList_Check( keywords ) && !PyTuple_Check( keywords ) )
        throw std::runtime_error( "Unexpected value for 'COMPORTEMENT'." );

    // Create syntax
    CommandSyntax cmdSt( "code_aster.Cata.Commons.c_comportement.C_COMPORTEMENT_SNL" );
    PyObject *kwfact = PyDict_New();
    PyDict_SetItemString( kwfact, "COMPORTEMENT", keywords );
    cmdSt.define( kwfact );

    // Build objects
    AS_ASSERT( _behavProp->build() );
    Py_DECREF( kwfact );
};
