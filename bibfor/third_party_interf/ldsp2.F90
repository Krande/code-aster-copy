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

subroutine ldsp2(pc, x1, y, ierr)
!
#include "asterf_petsc.h"
    use aster_petsc_module
    use petsc_data_module
!
    implicit none
! person_in_charge: natacha.bereux at edf.fr
!
#include "asterfort/amumph.h"
#include "asterfort/assert.h"

#ifdef ASTER_HAVE_PETSC
    complex(kind=8) :: cbid
    integer(kind=8) :: iret
    aster_logical :: prepos
!----------------------------------------------------------------
!     Variables PETSc
    PC, intent(in)              :: pc
    Vec, intent(in)             :: x1
    Vec, intent(inout)          :: y
    PetscErrorCode, intent(out) ::  ierr
!
    PetscScalar :: xx(1)
    PetscOffset :: xidx
!----------------------------------------------------------------
!
! --  COPIE DU VECTEUR D'ENTREE, CAR ERREUR S'IL EST TRANSFORME
    call VecCopy(x1, xlocal, ierr)
    ASSERT(ierr .eq. 0)
!
! --  RECUPERATION DES VALEURS DU VECTEUR SUR LES DIFFERENTS PROCS
    call VecScatterBegin(xscatt, xlocal, xglobal, INSERT_VALUES, SCATTER_FORWARD, ierr)
    ASSERT(ierr .eq. 0)
    call VecScatterEnd(xscatt, xlocal, xglobal, INSERT_VALUES, SCATTER_FORWARD, ierr)
    ASSERT(ierr .eq. 0)
!
    call VecGetArray(xglobal, xx, xidx, ierr)
    ASSERT(ierr .eq. 0)
!
! --  APPEL A LA ROUTINE DE PRECONDITIONNEMENT (DESCENTE/REMONTEE)
    cbid = dcmplx(0.d0, 0.d0)
    prepos = .true.
    call amumph('RESOUD', spsomu, spmat, xx(xidx+1), [cbid], &
                ' ', 1_8, iret, prepos)
!
! --  ENVOI DES VALEURS DU VECTEUR SUR LES DIFFERENTS PROCS
    call VecRestoreArray(xglobal, xx, xidx, ierr)
    ASSERT(ierr .eq. 0)
    call VecScatterBegin(xscatt, xglobal, y, INSERT_VALUES, SCATTER_REVERSE, ierr)
    ASSERT(ierr .eq. 0)
    call VecScatterEnd(xscatt, xglobal, y, INSERT_VALUES, SCATTER_REVERSE, ierr)
    ASSERT(ierr .eq. 0)
!
    ierr = to_petsc_int(iret)
!
#else
!
!     DECLARATION BIDON POUR ASSURER LA COMPILATION
    integer(kind=8) :: pc, x1, y, ierr
    integer(kind=8) :: idummy
    idummy = pc+x1+y
    ierr = 0
!
#endif
!
end subroutine
