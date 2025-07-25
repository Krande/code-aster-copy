! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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

subroutine elg_allocvr(vect1, n1)
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
    use aster_petsc_module
    implicit none
! person_in_charge: natacha.bereux@edf.fr
!
! BUT : allouer un vecteur PETSc réel de longueur n1
!---------------------------------------------------------------
#include "asterc/asmpi_comm.h"
#include "asterfort/utmess.h"
!
#ifdef ASTER_HAVE_PETSC
    Vec :: vect1
    integer(kind=8) :: n1
!
!================================================================
    PetscErrorCode :: ierr
    PetscInt :: bs
    mpi_int :: mpicomm
!----------------------------------------------------------------
    bs = 1
!    call asmpi_comm('GET_WORLD', mpicomm)
    mpicomm = PETSC_COMM_SELF
    call VecCreate(mpicomm, vect1, ierr)
    call VecSetBlockSize(vect1, bs, ierr)
    call VecSetType(vect1, VECSEQ, ierr)
    call VecSetSizes(vect1, PETSC_DECIDE, to_petsc_int(n1), ierr)
#else
    integer(kind=8) :: vect1, n1
    integer(kind=8) :: idummy
    call utmess('F', 'ELIMLAGR_1')
    idummy = vect1+n1
#endif
end subroutine
