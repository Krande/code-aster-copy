/**
 * @file LinearSolver.cxx
 * @brief Initialisation des renumeroteurs autorises pour les solvers
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "astercxx.h"

#include "Solvers/LinearSolver.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/ResultNaming.h"

LinearSolver::LinearSolver( const std::string name )
    : DataStructure( name, 19, "SOLVEUR" ), _isEmpty( true ),
      _charValues( JeveuxVectorChar24( getName() + ".SLVK" ) ),
      _doubleValues( JeveuxVectorReal( getName() + ".SLVR" ) ),
      _integerValues( JeveuxVectorLong( getName() + ".SLVI" ) ),
      _petscOptions( JeveuxVectorChar80( getName() + ".SLVO" ) ), _matrix( nullptr ),
      _matrixPrec(
          new AssemblyMatrixDisplacementReal( ResultNaming::getNewResultName() + ".PREC" ) ),
      _commandName( "SOLVEUR" ), _xfem( false ), _keywords( NULL ){

                                                 };

void LinearSolver::setKeywords( PyObject *user_keywords ) {
    _isEmpty = true;
    Py_XDECREF( _keywords );
    _keywords = user_keywords;
    Py_INCREF( _keywords );
#ifdef ASTER_DEBUG_CXX
    PYDBG( "setKeywords:", _keywords );
#endif
}

PyObject *LinearSolver::getKeywords() const {
    /* Returns a dict containing the SOLVEUR keywords.
     *
     * Return value: New reference.
     */
    AS_ASSERT( _keywords != NULL );
    PyObject *dict = PyDict_New();
    std::string mcf = "SOLVEUR";
    PyDict_SetItemString( dict, mcf.c_str(), _keywords );
    return dict;
}

bool LinearSolver::build() {
    if ( _charValues->exists() ) {
        _charValues->deallocate();
        _doubleValues->deallocate();
        _integerValues->deallocate();
    }
    std::string newName( getName() );
    newName.resize( 19, ' ' );

    // Definition du bout de fichier de commande pour SOLVEUR
    CommandSyntax cmdSt( _commandName );
    cmdSt.setResult( getName(), getType() );

    PyObject *dict = getKeywords();
    cmdSt.define( dict );

    std::string base( "G" );
    std::string xfem( "   " );
    if ( _xfem ) {
        xfem = "OUI";
    }
    CALLO_CRESOL_WRAP( newName, base, xfem );
    _isEmpty = false;
    Py_DECREF( dict );

    return true;
};

bool LinearSolver::deleteFactorizedMatrix() {

    if ( _matrix && _matrix->isFactorized() && get_sh_jeveux_status() == 1 ) {
        CALLO_DETMATRIX( _matrix->getName() );
        _matrix->deleteFactorizedMatrix();
    }

    return true;
};

bool LinearSolver::factorize( AssemblyMatrixDisplacementRealPtr currentMatrix ) {
    if ( _isEmpty )
        build();

    deleteFactorizedMatrix();
    _matrix = currentMatrix;

    const std::string solverName( getName() + "           " );
    std::string base( "G" );
    ASTERINTEGER cret = 0, npvneg = 0, istop = -9999;
    const std::string matpre( _matrixPrec->getName() );
    const std::string matass = _matrix->getName();

    // Definition du bout de fichier de commande pour SOLVEUR
    CommandSyntax cmdSt( _commandName );
    cmdSt.setResult( getName(), getType() );

    PyObject *dict = getKeywords();
    cmdSt.define( dict );

    CALLO_MATRIX_FACTOR( solverName, base, &cret, _matrixPrec->getName(), matass, &npvneg, &istop );

    _matrix->_isFactorized = true;
    _matrix->setSolverName( getSolverName() );

    Py_DECREF( dict );
    return true;
};

FieldOnNodesRealPtr LinearSolver::solve( const FieldOnNodesRealPtr currentRHS,
                                         FieldOnNodesRealPtr result ) const {

    if ( !_matrix ) {
        throw std::runtime_error( "Matrix must be factored first" );
    }

    if ( result->getName() == "" )
        result = boost::make_shared< FieldOnNodesReal >();

    try {
        if ( !result->getDOFNumbering() && currentRHS->getDOFNumbering() ) {
            result->setDOFNumbering( currentRHS->getDOFNumbering() );
        }
    } catch ( ... ) {
    }

    std::string blanc( " " );
    ASTERINTEGER nsecm = 0, istop = 0, iret = 0;
    ASTERDOUBLE rdummy = 0., cdummy = 0.;
    bool prepos( true );
    std::string base( JeveuxMemoryTypesNames[Permanent] );

    CALLO_RESOUD( _matrix->getName(), _matrixPrec->getName(), getName(), blanc, &nsecm,
                  currentRHS->getName(), result->getName(), base, &rdummy, &cdummy, blanc,
                  (ASTERLOGICAL *)&prepos, &istop, &iret );

    _matrix->setSolverName( getSolverName() );

    return result;
};

FieldOnNodesRealPtr LinearSolver::solveWithDirichletBC( const FieldOnNodesRealPtr currentRHS,
                                                        const FieldOnNodesRealPtr dirichletBCField,

                                                        FieldOnNodesRealPtr result ) const {

    if ( !_matrix ) {
        throw std::runtime_error( "Matrix must be factored first" );
    }

    if ( result->getName().empty() )
        result = boost::make_shared< FieldOnNodesReal >();

    try {
        if ( !result->getDOFNumbering() && currentRHS->getDOFNumbering() ) {
            result->setDOFNumbering( currentRHS->getDOFNumbering() );
        }
    } catch ( ... ) {
    }

    std::string blanc( " " );
    ASTERINTEGER nsecm = 0, istop = 0, iret = 0;
    ASTERDOUBLE rdummy = 0., cdummy = 0.;
    bool prepos( true );
    std::string base( JeveuxMemoryTypesNames[Permanent] );

    CALLO_RESOUD( _matrix->getName(), _matrixPrec->getName(), getName(),
                  dirichletBCField->getName(), &nsecm, currentRHS->getName(), result->getName(),
                  base, &rdummy, &cdummy, blanc, (ASTERLOGICAL *)&prepos, &istop, &iret );

    _matrix->setSolverName( getSolverName() );

    return result;
};
