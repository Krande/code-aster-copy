/**
 * @file PetscRedistribute.cxx
 * @brief Inspired from petsc PCTelescope, here is some wrapping functions
 *        that helps repartioning a matrix on a subcommunicator.
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
#include "aster_pybind.h"
#ifdef ASTER_HAVE_PETSC
#include "petsc.h"
#include "petscsystypes.h"
#endif

#ifdef ASTER_HAVE_PETSC
static PetscErrorCode redistribute_petsc( Mat mat, int subCommSize, Mat *new_mat );
#else
void redistribute_petsc();
#endif

py::object redistributePetscMat( py::object pMat, int subCommSize );
