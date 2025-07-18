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

module lmp_module
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
!
! person_in_charge: natacha.bereux at edf.fr
!
    use aster_petsc_module
    use lmp_context_type
!
    implicit none
!
    private
#include "asterf.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
!
    public :: lmp_apply_right, lmp_update, lmp_destroy
!
!
#ifdef ASTER_HAVE_PETSC
contains
!
! lmp_update est appelé APRES la résolution (KSPSolve)
!
    subroutine lmp_update(lmp_pc, ksp, ierr)
        !
        use lmp_data_module, only: lmp_context, lmp_is_setup, reac_lmp
        ! Dummy arguments
        PC, intent(inout)           :: lmp_pc
        KSP, intent(in)             :: ksp
        PetscErrorCode, intent(out) :: ierr
        ! Local variables
        PetscInt :: maxits
        integer(kind=8) :: ifm, niv
        aster_logical :: verbose
        !
        call infniv(ifm, niv)
        verbose = niv == 2

        call KSPGetIterationNumber(ksp, maxits, ierr)
        ASSERT(ierr == 0)
        if (verbose) then
            call utmess('I', 'PETSC_24', si=to_aster_int(maxits))
        end if
        if (lmp_is_setup) then
            if (maxits > reac_lmp) then
                ! todo : il faudrait normalement détruire le LMP, mais cela fonctionne moins bien
                ! à comprendre ...
                !    call free_lmp_context( lmp_context )
                !    lmp_is_setup =.false.
            end if
            if (verbose) then
                call utmess('I', 'PETSC_27', si=reac_lmp)
            end if
        end if
        if (.not. lmp_is_setup) then
            if (maxits > reac_lmp) then
                ! On construit un nouveau lmp
                if (verbose) then
                    call utmess('I', 'PETSC_26')
                end if
                lmp_context = new_lmp_context()
                call build_lmp_context(ksp, lmp_context)
                lmp_is_setup = .true.
            end if
        end if
        !
    end subroutine lmp_update
!
! PCApply (Right) operation used for the PCShell based on the LMP Method
!----------------------------------------------------------------
!
!   APPLICATION DU LMP: y=(I+zz_1*yy_1^T)...(I+zz_lmp*yy_lmp^T) x
!
! ---------------------------------------------------------------
!
    subroutine lmp_apply_right(lmp_pc, x, y, ierr)
        !
        use lmp_data_module, only: lmp_context, lmp_is_setup
        !
        ! Dummy arguments
        !
        PC, intent(in)    :: lmp_pc
        Vec, intent(in)   :: x
        Vec, intent(inout):: y
        PetscErrorCode, intent(out) :: ierr
        !
        ! Local Variables
        PetscInt :: jj
        !
        ! --  Copie du vecteur d'entree
        call VecCopy(x, y, ierr)
        if (lmp_is_setup) then
            ! --  sk=yy^T x
            do jj = 1, lmp_context%ritzeff
                call VecDot(lmp_context%yy(jj), y, lmp_context%sk(jj), ierr)
                ASSERT(ierr == 0)
            end do
! --  y=I+zz*sk
            call VecMAXPY(y, lmp_context%ritzeff, lmp_context%sk, lmp_context%zz, ierr)
            ASSERT(ierr == 0)
        end if
        !
    end subroutine lmp_apply_right
!
    subroutine lmp_destroy(lmp_pc, ierr)
        !
        use lmp_data_module, only: lmp_context, lmp_is_setup
        !
        ! Dummy arguments
        !
        PC, intent(inout)    :: lmp_pc
        PetscErrorCode, intent(out) :: ierr
        !
        call free_lmp_context(lmp_context)
        lmp_is_setup = .false.
        !
        ierr = 0
        !
    end subroutine lmp_destroy
!
#else
contains
!
    subroutine lmp_update(lmp_pc, ksp, ierr)
        !
        ! Dummy arguments
        integer(kind=8)  :: lmp_pc
        integer(kind=8)  :: ksp
        integer(kind=8)  :: ierr
        lmp_pc = 0
        ksp = 0
        ierr = 0
        ASSERT(.false.)
    end subroutine lmp_update
!
    subroutine lmp_apply_right(lmp_pc, x, y, ierr)
        integer(kind=8) :: lmp_pc
        integer(kind=8) :: x
        integer(kind=8) :: y
        integer(kind=8) :: ierr
        lmp_pc = 0
        x = 0
        y = 0
        ierr = 0
        ASSERT(.false.)
    end subroutine lmp_apply_right
!
    subroutine lmp_destroy(lmp_pc, ierr)
        integer(kind=8) :: lmp_pc
        integer(kind=8) :: ierr
        lmp_pc = 0
        ierr = 0
        ASSERT(.false.)
    end subroutine lmp_destroy
#endif
end module lmp_module
