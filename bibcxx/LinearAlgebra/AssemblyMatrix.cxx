/**
 * @file AssemblyMatrix.cxx
 * @brief Implementation de AssemblyMatrix vide car AssemblyMatrix est un template
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

#include "LinearAlgebra/AssemblyMatrix.h"

#include "Solvers/LinearSolver.h"

// Specialization for <double, Displacement>
template <>
void AssemblyMatrix< ASTERDOUBLE, Displacement >::setValues( const VectorLong &idx,
                                                             const VectorLong &jdx,
                                                             const VectorReal &values ) {
    const ASTERINTEGER dim = idx.size();
    if ( idx.size() != jdx.size() || idx.size() != values.size() ) {
        throw std::runtime_error( "All lists must have same length" );
    }
    CALLO_MATR_ASSE_SET_VALUES( getName(), &dim, idx.data(), jdx.data(), values.data() );
    _isFactorized = false;
    _matrixValues->build( true );
};

template <>
void AssemblyMatrix< ASTERDOUBLE, Displacement >::scale( const VectorReal &lvect,
                                                         const VectorReal &rvect ) {
    const ASTERINTEGER matSize = size( true )[0];
    if ( lvect.size() != matSize || rvect.size() != matSize ) {
        throw std::runtime_error( "Arguments must share the matrix size" );
    }
    CALLO_MATR_ASSE_SCALE( getName(), lvect.data(), rvect.data() );
    _isFactorized = false;
};

template <>
void AssemblyMatrix< ASTERDOUBLE, Displacement >::defineSolver() {
    _solver = std::make_shared< LinearSolver >( ljust( getName(), 8 ) + ".SOLVEUR   " );
}

template <>
void AssemblyMatrix< ASTERDOUBLE, Displacement >::applyDirichletBC(
    const FieldOnNodesReal &DirichletBC, FieldOnNodesReal &Rhs ) const {
    if ( get_sh_jeveux_status() == 1 ) {
        MATR_ASSE_COMPUTE_KINEMATIC_RHS( getName(), DirichletBC.getName(), Rhs.getName() );
    } else {
        raiseAsterError( "Jeveux is not ready!" );
    }
};

template <>
void AssemblyMatrix< ASTERCOMPLEX, Displacement >::defineSolver() {
    _solver = std::make_shared< LinearSolver >( ljust( getName(), 8 ) + ".SOLVEUR   " );
}

// Specialization for <double, Temperature>
template <>
void AssemblyMatrix< ASTERDOUBLE, Temperature >::setValues( const VectorLong &idx,
                                                            const VectorLong &jdx,
                                                            const VectorReal &values ) {
    const ASTERINTEGER dim = idx.size();
    if ( idx.size() != jdx.size() || idx.size() != values.size() ) {
        throw std::runtime_error( "All lists must have same length" );
    }
    CALLO_MATR_ASSE_SET_VALUES( getName(), &dim, idx.data(), jdx.data(), values.data() );
    _isFactorized = false;
    _matrixValues->build( true );
};

template <>
void AssemblyMatrix< ASTERDOUBLE, Temperature >::scale( const VectorReal &lvect,
                                                        const VectorReal &rvect ) {
    const ASTERINTEGER matSize = size( true )[0];
    if ( lvect.size() != matSize || rvect.size() != matSize ) {
        throw std::runtime_error( "Arguments must share the matrix size" );
    }
    CALLO_MATR_ASSE_SCALE( getName(), lvect.data(), rvect.data() );
    _isFactorized = false;
};

template <>
void AssemblyMatrix< ASTERDOUBLE, Temperature >::defineSolver() {
    _solver = std::make_shared< LinearSolver >( ljust( getName(), 8 ) + ".SOLVEUR   " );
}

template <>
void AssemblyMatrix< ASTERDOUBLE, Temperature >::applyDirichletBC(
    const FieldOnNodesReal &DirichletBC, FieldOnNodesReal &Rhs ) const {
    if ( get_sh_jeveux_status() == 1 ) {
        MATR_ASSE_COMPUTE_KINEMATIC_RHS( getName(), DirichletBC.getName(), Rhs.getName() );
    } else {
        raiseAsterError( "Jeveux is not ready!" );
    }
};

// Specialization for <double, Pressure>
template <>
void AssemblyMatrix< ASTERDOUBLE, Pressure >::setValues( const VectorLong &idx,
                                                         const VectorLong &jdx,
                                                         const VectorReal &values ) {
    const ASTERINTEGER dim = idx.size();
    if ( idx.size() != jdx.size() || idx.size() != values.size() ) {
        throw std::runtime_error( "All lists must have same length" );
    }
    CALLO_MATR_ASSE_SET_VALUES( getName(), &dim, idx.data(), jdx.data(), values.data() );
    _isFactorized = false;
    _matrixValues->build( true );
};

template <>
void AssemblyMatrix< ASTERDOUBLE, Pressure >::applyDirichletBC( const FieldOnNodesReal &DirichletBC,
                                                                FieldOnNodesReal &Rhs ) const {
    if ( get_sh_jeveux_status() == 1 ) {
        MATR_ASSE_COMPUTE_KINEMATIC_RHS( getName(), DirichletBC.getName(), Rhs.getName() );
    } else {
        raiseAsterError( "Jeveux is not ready!" );
    }
};

template <>
void AssemblyMatrix< ASTERCOMPLEX, Pressure >::defineSolver() {
    _solver = std::make_shared< LinearSolver >( ljust( getName(), 8 ) + ".SOLVEUR   " );
}
