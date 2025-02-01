/**
 * @file LinearSolver.cxx
 * @brief Initialisation des renumeroteurs autorises pour les solvers
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

#include "aster_pybind.h"

#include "Supervis/CommandSyntax.h"
#include "Supervis/ResultNaming.h"

LinearSolver::LinearSolver( const std::string name )
    : DataStructure( name, 19, "SOLVEUR" ),
      _isBuilt( false ),
      _charValues( JeveuxVectorChar24( getName() + ".SLVK" ) ),
      _doubleValues( JeveuxVectorReal( getName() + ".SLVR" ) ),
      _integerValues( JeveuxVectorLong( getName() + ".SLVI" ) ),
      _petscOptions( JeveuxVectorChar80( getName() + ".SLVO" ) ),
      _matrix( nullptr ),
      _matrixPrec( nullptr ),
      _cataPath( std::string() ),
      _xfem( false ),
      _keywords( py::none() ) {

      };

void LinearSolver::setKeywords( py::object &user_keywords ) {
    _isBuilt = false;
    _keywords = user_keywords;
}

py::dict LinearSolver::getKeywords() const {
    /* Returns a dict containing the SOLVEUR keywords.
     *
     * Return value: New reference.
     */
    AS_ASSERT( !_keywords.is_none() );
    py::dict dict( py::arg( "SOLVEUR" ) = _keywords );
    return dict;
}

bool LinearSolver::build() {
    if ( _charValues.exists() ) {
        _charValues->deallocate();
        _doubleValues->deallocate();
        _integerValues->deallocate();
        _petscOptions->deallocate();
    }
    std::string newName( getName() );
    newName.resize( 19, ' ' );

    // Definition du bout de fichier de commande pour SOLVEUR
    CommandSyntax cmdSt( _cataPath );
    cmdSt.setResult( getName(), getType() );

    cmdSt.define( getKeywords() );

    std::string base( "G" );
    std::string xfem( "   " );
    if ( _xfem ) {
        xfem = "OUI";
    }
    CALLO_CRESOL_WRAP( newName, base, xfem );
    _isBuilt = true;

    return true;
};

bool LinearSolver::deleteFactorizedMatrix() {

    if ( _matrix && _matrix->isFactorized() && get_sh_jeveux_status() == 1 ) {
        _matrix->deleteFactorizedMatrix();
        CALLO_DETMATRIX( _matrix->getName() );
    }

    return true;
};

bool LinearSolver::factorize( const BaseAssemblyMatrixPtr currentMatrix, bool raiseException ) {

    deleteFactorizedMatrix();
    _matrix = currentMatrix;

    if ( _matrix->getMesh()->isParallel() ) {
        if ( !supportParallelMesh() ) {
            raiseAsterError( "Solver does not support a parallel mesh" );
        }
    }

    if ( !_isBuilt )
        build();

    const std::string solverName( getName() + "           " );
    std::string base( "G" );
    ASTERINTEGER cret = 0, npvneg = 0;
    ASTERINTEGER istop = -9999;
    if ( raiseException ) {
        istop = 2;
    }

    _matrixPrec = _matrix->getEmptyMatrix( ResultNaming::getNewResultName() );
    const std::string matpre( _matrixPrec->getName() );
    const std::string matass( _matrix->getName() );

    // Definition du bout de fichier de commande pour SOLVEUR
    CommandSyntax cmdSt( _cataPath );
    cmdSt.setResult( getName(), getType() );

    cmdSt.define( getKeywords() );

    CALLO_MATRIX_FACTOR( solverName, base, &cret, _matrixPrec->getName(), matass, &npvneg, &istop );
    if ( raiseException && cret != 0 ) {
        throw SolverErrorCpp( "FACTOR_13" );
    }

    _matrix->setFactorized( true );
    _matrix->setSolverName( getSolverName() );
    if ( getSolverName() == "GCPC" )
        _matrixPrec->updateDOFNumbering();

    return true;
};

void LinearSolver::_solve( const std::string &rhsName, const std::string &diriName,
                           const std::string &resultName ) const {
    std::string blanc( " " );
    ASTERINTEGER nsecm = 0, istop = 0, iret = 0;
    ASTERDOUBLE rdummy = 0., cdummy = 0.;
    bool prepos( true );
    std::string base( JeveuxMemoryTypesNames[Permanent] );

    CALLO_RESOUD( _matrix->getName(), _matrixPrec->getName(), getName(), diriName, &nsecm, rhsName,
                  resultName, base, &rdummy, &cdummy, blanc, (ASTERLOGICAL *)&prepos, &istop,
                  &iret );
}

FieldOnNodesRealPtr LinearSolver::solve( const FieldOnNodesRealPtr currentRHS,
                                         const FieldOnNodesRealPtr dirichletBCField ) const {

    if ( !_matrix ) {
        raiseAsterError( "Matrix must be factored first" );
    }

    auto result = std::make_shared< FieldOnNodesReal >( _matrix->getDOFNumbering() );

    std::string diriName( " " );
    if ( dirichletBCField )
        diriName = dirichletBCField->getName();

    _solve( currentRHS->getName(), diriName, result->getName() );

    return result;
};

FieldOnNodesComplexPtr LinearSolver::solve( const FieldOnNodesComplexPtr currentRHS,
                                            const FieldOnNodesComplexPtr dirichletBCField ) const {

    if ( !_matrix ) {
        raiseAsterError( "Matrix must be factored first" );
    }

    auto result = std::make_shared< FieldOnNodesComplex >( _matrix->getDOFNumbering() );

    std::string diriName( " " );
    if ( dirichletBCField )
        diriName = dirichletBCField->getName();

    _solve( currentRHS->getName(), diriName, result->getName() );

    return result;
};
