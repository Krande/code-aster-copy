/**
 * @file PhysicalProblem.cxx
 * @brief Implementation of class PhysicalProblem
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

#include "Studies/PhysicalProblem.h"

#include "aster_pybind.h"

#include "Numbering/DOFNumbering.h"
#include "Numbering/ParallelDOFNumbering.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/Exceptions.h"

PhysicalProblem::PhysicalProblem( const ModelPtr curModel, const MaterialFieldPtr curMat,
                                  const ElementaryCharacteristicsPtr cara )
    : _model( curModel ),
      _mesh( curModel->getMesh() ),
      _materialField( curMat ),
      _elemChara( cara ),
      _listOfLoads( std::make_shared< ListOfLoads >( _model ) ),
      _dofNume( nullptr ),
      _codedMater( nullptr ),
      _behavProp( nullptr ),
      _externVarRefe( nullptr ) {

    // Add checks
    if ( _elemChara ) {
        if ( _model != _elemChara->getModel() ) {
            const std::string msg = "Inconsistent models: " + _model->getName() + " vs " +
                                    _elemChara->getModel()->getName();
            AS_ABORT( msg );
        }
    }

    if ( _materialField ) {
        if ( _mesh != _materialField->getMesh() ) {
            const std::string msg = "Inconsistent meshes: " + _mesh->getName() + " vs " +
                                    _materialField->getMesh()->getName();
            AS_ABORT( msg );
        }

        _codedMater = std::make_shared< CodedMaterial >( _materialField, _model );
        _codedMater->allocate( true );
    }
};

void PhysicalProblem::setDOFNumbering( const BaseDOFNumberingPtr dofNume ) {
    if ( dofNume->getMesh() != _mesh ) {
        const std::string msg =
            "Inconsistent mesh: " + _mesh->getName() + " vs " + dofNume->getMesh()->getName();
        AS_ABORT( msg );
    }

    auto model = dofNume->getModel();
    if ( model && model != _model ) {
        const std::string msg =
            "Inconsistent models: " + _model->getName() + " vs " + model->getName();
        AS_ABORT( msg );
    }

    auto listOfLoads = dofNume->getListOfLoads();
    if ( listOfLoads && !listOfLoads->isEmpty() && listOfLoads != _listOfLoads ) {
        const std::string msg = "Inconsistent list of loads: " + _listOfLoads->getName() + " vs " +
                                listOfLoads->getName();
        AS_ABORT( msg );
    }

    _dofNume = dofNume;
};

bool PhysicalProblem::computeDOFNumbering() {
    // create dofNume
#ifdef ASTER_HAVE_MPI
    if ( _mesh->isParallel() )
        _dofNume = std::make_shared< ParallelDOFNumbering >();
    else
#endif /* ASTER_HAVE_MPI */
        _dofNume = std::make_shared< DOFNumbering >();

    _dofNume->setModel( _model );
    _dofNume->setListOfLoads( _listOfLoads );

    return _dofNume->computeNumbering();
};

void PhysicalProblem::computeBehaviourProperty( py::object &keywords,
                                                const std::string &initialState,
                                                const ASTERINTEGER verbosity ) {
    // Create object for behaviour
    _behavProp = std::make_shared< BehaviourProperty >( _model, _materialField );
    _behavProp->setInitialState( initialState == "OUI" );
    _behavProp->setVerbosity( verbosity > 1 );

    // Check input PyObject
    if ( !PyDict_Check( keywords.ptr() ) && !PyList_Check( keywords.ptr() ) &&
         !PyTuple_Check( keywords.ptr() ) )
        throw std::runtime_error( "Unexpected value for 'COMPORTEMENT'." );

    // Create syntax
    CommandSyntax cmdSt( "code_aster.Cata.Commons.c_comportement.C_COMPORTEMENT_SNL" );
    py::dict kwfact( py::arg( "COMPORTEMENT" ) = keywords );
    cmdSt.define( kwfact );

    // Build objects
    AS_ASSERT( _behavProp->build() );
};

FieldOnCellsRealPtr PhysicalProblem::getExternalStateVariables( const ASTERDOUBLE &time ) const {

    // Create field
    auto FEDesc = getModel()->getFiniteElementDescriptor();
    auto field = std::make_shared< FieldOnCellsReal >( FEDesc );

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( getModel()->getName(), 24 );
    std::string materialFieldName = ljust( getMaterialField()->getName(), 24 );
    auto currElemChara = getElementaryCharacteristics();
    std::string elemCharaName( " " );
    if ( currElemChara )
        elemCharaName = currElemChara->getName();
    elemCharaName.resize( 24, ' ' );
    std::string fieldName = ljust( field->getName(), 19 );

    // Output
    std::string out( ' ', 2 );
    std::string base( "G" );

    // Call Fortran WRAPPER
    CALLO_VRCINS_WRAP( modelName, materialFieldName, elemCharaName, &time, fieldName, out, base );

    return field;
}

void PhysicalProblem::computeReferenceExternalStateVariables() {

    // Create field
    auto FEDesc = getModel()->getFiniteElementDescriptor();
    _externVarRefe = std::make_shared< FieldOnCellsReal >( FEDesc );

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( getModel()->getName(), 8 );
    std::string materialFieldName = ljust( getMaterialField()->getName(), 8 );
    auto currElemChara = getElementaryCharacteristics();
    std::string elemCharaName( " ", 8 );
    if ( currElemChara )
        elemCharaName = std::string( currElemChara->getName(), 0, 8 );
    std::string fieldName = ljust( _externVarRefe->getName(), 19 );
    std::string base( "G" );

    // Call Fortran WRAPPER
    CALLO_VRCREF( modelName, materialFieldName, elemCharaName, fieldName, base );
}