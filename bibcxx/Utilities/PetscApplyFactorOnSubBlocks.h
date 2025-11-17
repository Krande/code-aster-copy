/**
 * @file PetscRedistribute.cxx
 * @brief Given a factored matrix, here is a petsc wrapping function
 *        to effciently compute from a sparse RHS matrixonly a subset
 *        of entries of solution matrix.
 * @author Nicolas Tardieu
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
#pragma once
#include "aster_pybind.h"
#ifdef ASTER_HAVE_PETSC
#include "petsc.h"
#include "petscsystypes.h"
#endif

#ifdef ASTER_HAVE_PETSC
PetscErrorCode MatSparseSolve_petsc( Mat FctMat, Mat RhsMat, IS ISet, Mat *SolMat,
                                     IS *JSet = NULL );
#else
void MatSparseSolve_petsc();
#endif

py::object applyFactorOnSubBlocks( py::object pyFctMat, py::object pyRhsMat, py::object pyISet,
                                   py::object pyJSet = py::none() );
