/**
 * @file PetscRedistributeInterface.cxx
 * @brief Interface python de PetscRedistribute
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

#include "PythonBindings/PetscRedistributeInterface.h"

#include "Utilities/PetscRedistribute.h"

void exportPetscRedistributeToPython( py::module_ &mod ) {
    mod.def( "redistributePetscMat", &redistributePetscMat,
             R"(Given a distributed matrix on an MPI communicator,
     this function returns a redistributed matrix on a subcommunicator.

     Arguments:
        pMat: the petsc4py matrix to redistribute.
        subCommSize: the size of the subcommunicator
     Outputs:
        outMat: the redistributed petsc4py matrix on a subcommunicator of size
                 subCommSize.
   )",
             py::arg( "pMat" ), py::arg( "subCommSize" ) );
}
