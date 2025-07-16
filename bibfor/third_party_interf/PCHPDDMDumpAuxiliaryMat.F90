! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF RD - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------

subroutine PCHPDDMDumpAuxiliaryMat(pc, is, aux)
    !
#include "asterf_types.h"
#include "asterf_petsc.h"
    !
    !
    use aster_petsc_module
    implicit none
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
    !----------------------------------------------------------------
    !
    !  CONSTRUCTION DE LA MATRICE DE RIGIDITE LOCALE DU SOUS-DOMAINE
    !
    !  En entrée : la matrice ASTER complète
    !  En sortie : La matrice PETSc définie dans l'appelant
    !
    !  Rq :
    !  - la matrice construite comporte tous les degrés de liberté du sous-domaine
    !    y compris les ghosts
    !  - dans la terminologie de la decomposition de domaine, la matrice construite
    !    est la matrice de Neumann
    !----------------------------------------------------------------
    !
#ifdef ASTER_HAVE_PETSC
    !
    !----------------------------------------------------------------
    !     Variables PETSc
    Mat            :: aux, P
    PC             :: pc
    IS             :: iss, is
    PetscViewer    :: viewer
    PetscInt       :: sizes(4)
    PetscMPIInt    :: rank, size
    character(len=100) ::  dir, name
    PetscErrorCode :: ierr
    PetscBool      :: set
    !----------------------------------------------------------------
    call asmpi_info(comm=PETSC_COMM_WORLD, rank=rank, size=size)
    dir = "/tmp"
   call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dump_dir", dir, set, ierr)
    ASSERT(ierr .eq. 0)
    write (name, "(A,A,I0.3,A,I0.3,A)") trim(dir), "/sizes_", rank, "_", size, ".dat"
    call PetscViewerBinaryOpen(PETSC_COMM_SELF, name, FILE_MODE_WRITE, viewer, ierr)
    ASSERT(ierr .eq. 0)
    call PCGetOperators(pc, PETSC_NULL_MAT, P, ierr)
    ASSERT(ierr .eq. 0)
    call MatGetLocalSize(P, sizes(1), sizes(2), ierr)
    ASSERT(ierr .eq. 0)
    call MatGetSize(P, sizes(3), sizes(4), ierr)
    ASSERT(ierr .eq. 0)
    call ISCreateGeneral(PETSC_COMM_SELF, to_petsc_int(4), sizes, PETSC_USE_POINTER, iss, ierr)
    ASSERT(ierr .eq. 0)
    call ISView(iss, viewer, ierr)
    ASSERT(ierr .eq. 0)
    call PetscViewerDestroy(viewer, ierr)
    ASSERT(ierr .eq. 0)
    call ISDestroy(iss, ierr)
    ASSERT(ierr .eq. 0)
    write (name, '(A,A,I0.3,A)') trim(dir), "/P_", size, ".dat"
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, name, FILE_MODE_WRITE, viewer, ierr)
    ASSERT(ierr .eq. 0)
    call MatView(P, viewer, ierr)
    ASSERT(ierr .eq. 0)
    call PetscViewerDestroy(viewer, ierr)
    ASSERT(ierr .eq. 0)
    write (name, '(A,A,I0.3,A,I0.3,A)') trim(dir), "/Neumann_", rank, "_", size, ".dat"
    call PetscViewerBinaryOpen(PETSC_COMM_SELF, name, FILE_MODE_WRITE, viewer, ierr)
    ASSERT(ierr .eq. 0)
    call MatView(aux, viewer, ierr)
    ASSERT(ierr .eq. 0)
    call PetscViewerDestroy(viewer, ierr)
    ASSERT(ierr .eq. 0)
    write (name, '(A,A,I0.3,A,I0.3,A)') trim(dir), "/is_", rank, "_", size, ".dat"
    call PetscViewerBinaryOpen(PETSC_COMM_SELF, name, FILE_MODE_WRITE, viewer, ierr)
    ASSERT(ierr .eq. 0)
    call ISView(is, viewer, ierr)
    ASSERT(ierr .eq. 0)
    call PetscViewerDestroy(viewer, ierr)
    ASSERT(ierr .eq. 0)
#else
    integer(kind=8) :: idummy
    integer(kind=8) :: pc, is, aux
#endif
    !
end subroutine
