/**
 * @file PetscApplyFactorOnSubBlocksInterface.cxx
 * @brief Interface python de applyFactorOnSubBlocks
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
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

#include "PythonBindings/PetscApplyFactorOnSubBlocksInterface.h"

#include "Utilities/PetscApplyFactorOnSubBlocks.h"

void exportApplyFactorOnSubBlocksToPython( py::module_ &mod ) {
    mod.def( "applyFactorOnSubBlocks", &applyFactorOnSubBlocks,
             R"(Given a MUMPS factor matrix and a sparse PETSc Mat as the right-hand side,
            this routine solves for all columns. If an index set I is provided, it first
            determines the set J of row indices where the RHS has nonzero entries, and
            returns a sparse MATAIJ solution containing only the IxJ block. Supplying I
            can speed up the solve, but is most effective when the RHS column count is
            comparable to or larger than |J|.
    Arguments:
        pyFctMat (PETSc.Mat): the MUMPS factor matrix.
        pyRHS (PETSc.Mat): the sparse MATSEQAIJ holding the multiple right hand sides
        pyISet (PETSc.IS): the index set for entries to compute
        pyJSet (PETSc.IS): optional off diagonal layout index set
    Outputs:
        outMat: the petsc4py aij matrix matrix holding the solutions
)",
             py::arg( "pyFctMat" ), py::arg( "pyRHS" ), py::arg( "pyISet" ),
             py::arg( "pyJSet" ) = py::none() );
}
